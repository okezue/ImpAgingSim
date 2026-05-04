> [!abstract] Thesis
> LLM inference currently fetches, dequantizes, or computes all bits of all model state by default, which is wasteful. Instead, it should keep a hardware-native low-bit prefix resident, compute certified uncertainty bounds from that prefix, and fetch refinement bits only for the pages, experts, or tokens whose uncertainty can still affect the next token.

> [!example] Connection to `Serif` and the dequantization wall
`Serif` shows that INT4 KV compression should reduce memory traffic, but on fixed-datapath accelerators the unpack/extend/scale path erases the benefit, with measured GPU overhead up to 409%, Trainium INT4 achieving only 53–57% of BF16 throughput despite smaller payloads, and FPGA DP4A fusion giving a 7.6× same-bitstream speedup for a 64×64 INT4 GEMM. TurboQuant and RaBitQ are strong baselines because they give principled vector quantization and inner-product/distance guarantees, but both are fundamentally about quantizing vectors or estimating distances, while PRISM focuses on **query-adaptive, certifiable, progressively refined inference state**.

---
## Core contributions
1. [[#Define certified progressive inference as a new problem]]
2. [[#PRISM-KV, a certified progressive KV cache]]
3. [[#PRISM-Gate, a certified progressive MoE router]]
4. [[#Native-prefix GPU implementation]]
5. [[#Attention softmax-mass certificates]]
6. [[#Output-error certificates]]
7. [[#Attention-read-aware progressive KV encoding]]
8. [[#Read-aware bit allocation]]
9. [[#Heterogeneous-memory scheduler]]
10. [[#End-to-end demo]]
11. [[#Negative result against static quantization]]
12. [[#Bridge to `Serif` and custom hardware]]
---
### Define certified progressive inference as a new problem
I introduce a new abstraction where each KV page, MoE gate block, expert-weight tile, or activation tile is stored as a nested bitstream $c^{(0)}(x), c^{(1)}(x), \ldots, c^{(M)}(x),$ where stage 0 is hardware-native and resident in fast memory, while later stages are refinements stored in slower memory. After reading stages $(0,\ldots,m)$, the runtime reconstructs $\hat{x}^{(m)}$ and knows a certified residual radius $|x-\hat{x}^{(m)}|_2 \le \rho_m.$
Here, our inference system asks a fundamental question: **is the current prefix enough to certify the next computation?** 
If yes, it stops. 
If no, it fetches another refinement stage only for the uncertain blocks. 
Basically, we treat quantization as a runtime decision protocol rather than a static procedure. 
### PRISM-KV, a certified progressive KV cache
Store each KV page as a native low-bit prefix plus residual bitplanes. The attention kernel first computes coarse score intervals from the prefix, then refines only pages whose upper-bound softmax mass remains non-negligible, and then, instead of reading the entire KV cache at high precision every token, the model reads a small always-resident prefix and only occasionally touches deeper memory.
### PRISM-Gate, a certified progressive MoE router
Keep coarse routing/gating information for all experts in fast memory. Certify the top-$k$ experts from intervals whenever possible and then fetch or refine only ambiguous experts near the routing boundary. This is potentially bigger than KV compression, because MoE models often have huge total parameter sets but activate only a few experts per token. 
### Native-prefix GPU implementation
For GPUs, the stage-0 code should use whatever the target hardware consumes efficiently without full scalar dequantization: FP8/INT8 on older targets, FP4-like block-scaled formats where available, and INT4/INT2 variants for custom backends. The crucial point is that the prefix should be **computable**and so refinements can be applied as residual corrections, but the always-read path must avoid the dequant wall from the `Serif` study.
### Attention softmax-mass certificates
For every page or token, PRISM maintains a lower and upper bound on its attention score and from these score intervals, it derives an upper bound on the total possible softmax mass of all ignored/unrefined pages, which produces a rigorous stopping rule where if the omitted mass is below tolerance, do not fetch more bits.
### Output-error certificates
Instead of stopping at a key-score error we can certify the actual attention output error via certificate combining omitted softmax mass, approximate-score perturbation among retained pages, and value reconstruction error.
### Attention-read-aware progressive KV encoding
Introduce a write-time encoder for PRISM that exploits the write-once/read-many structure of the KV cache. Each KV page is stored as a hardware-native low-bit prefix plus residual refinement stages. The prefix is not produced by ordinary per-token MSE quantization; it is optimized for future attention reads. For values, we derive a causal residual-feedback encoder whose readout error telescopes and is bounded by the total variation of the future attention weights. For keys, we derive a first-order softmax-sensitivity objective showing that key precision should be allocated according to $p_{\tau,t}(v_t-o_\tau)q_{\tau,j}$, not merely according to key norm or quantization MSE. This encoder makes PRISM’s resident prefix more informative, reducing the frequency with which runtime certificates must fetch deeper refinement bits.

>[!info] PRISM's write-time encoder
>The purpose of this is to make the resident prefix useful enough that PRISM rarely fetches refinement bits.
>The static encoder should produce: $$x = \hat{x}^{(0)} + r^{(1)} + r^{(2)} + \cdots + r^{(M)},$$  where $\hat{x}^{(0)}$ is the always-resident low-bit prefix, and $r^{(m)}$ are progressively stored residual stages. 
>After stage $m$, PRISM knows: $$\left|x-\hat{x}^{(m)}\right|_2 \le \rho_m,$$ and the runtime then uses $\rho_m$ to decide whether more precision is needed.
>For values, suppose a scalar value stream along token order is $x_1,\ldots,x_T$. Instead of independently rounding each $x_t$, the write-time encoder uses residual feedback: $$u_t = x_t + r_{t-1},$$ $$\hat{x}_t = Q_b(u_t),$$ $$r_t = u_t-\hat{x}_t.$$ Then the quantization error is $\varepsilon_t = \hat{x}*t-x_t = r*{t-1}-r_t.$ 
>Now consider a future attention read $y = \sum_{t=1}^{T} a_t x_t$, where the readout error is $$\hat{y}-y=\sum_{t=1}^{T} a_t(r_{t-1}-r_t),$$ which telescopes: $$\hat{y}-y=\sum_{t=1}^{T-1}(a_{t+1}-a_t)r_t-a_T r_T,$$ so if $|r_t|\le \Delta/2$, $$|\hat{y}-y|\le\frac{\Delta}{2}\left(\operatorname{TV}(a)+|a_T|\right),$$ where $\operatorname{TV}(a)=\sum_{t=1}^{T-1}|a_{t+1}-a_t|.$ 
>Then for value-cache attention, $a_t=p_{\tau,t},$ the future attention weights, so $$|\delta o_{\tau,c}^{(V)}|\le\frac{\Delta_c}{2}\left(\operatorname{TV}(p_\tau)+p_{\tau,T}\right),$$ which is our theorem that **the value-prefix error is controlled by how much the future attention row varies along token order, not by raw MSE alone.** This is why the resident low-bit prefix can be much more accurate than an ordinary low-bit quantizer when attention is smooth or distributed.
>#### For keys, it's softmax-senstivity-aware coding
>Let $$s_{\tau,t}=\frac{q_\tau^\top k_t}{\sqrt d}, \qquad p_tau = \operatorname{softmax}(s_\tau),\qquad o_\tau=\sum_t p_{\tau, t}v_t,$$ then for small key/value perturbations, $$\delta o_\tau \approx \sum_t p_{\tau,t}\delta v_t + \sum_t p_{\tau,t}(v_t-o_\tau)\frac{q_\tau^\top \delta} k_t{\sqrt d},$$ which means value error is weighted by $p_{\tau,t}$ why key error is weighted by $p_{\tau,t}(v_t-o_\tau)q_{\tau,j}$, so the write-time encoder should allocate more prefix precision to key coordinates and pages where this first-order sensitivity is large, and fewer bits where key perturbations have little effect on the output.
### Read-aware bit allocation
Instead of allocating bits based on tensor MSE, do it according to certified downstream sensitivity so then heads/pages with low attention variation and large margins get fewer refinement bits and sink heads, retrieval heads, and uncertain MoE routers get more. This gives a theoretically derived mixed-precision policy.
### Heterogeneous-memory scheduler
A scheduler that places stage 0 in GPU HBM, stage 1 in slower GPU memory or host-pinned memory, stage 2 in CPU/CXL memory, and cold residuals deeper still and helps minimize expected bytes moved per token subject to the certificate.
### End-to-end demo
Demo running a context length or model size that does not fit in normal HBM at useful throughput by keeping only coarse prefixes resident reporting accuracy, largest mode/context served on fixed hardware, bytes fetched per token from each tier, inter-token latency, and throughput.
### Negative result against static quantization
A core proof point is that static 2-bit/3-bit/4-bit quantization wastes bits on inactive state, as static quantization pays for all blocks every token, while PRISM pays for uncertainty only.
### Bridge to `Serif` and custom hardware
Multiple physical instantiations for PRISM: on commodity GPUs, it avoids reading/fetching deeper bits unless needed, on FPGA/ASIC-style backends, `Serif`-style fused sub-byte units make the prefix path truly dequanti-free, which lets the same algorithmic abstraction span diversified hardware. Hopefully `Serif` results would have already been published at this point.

---

> [!important]
> ## Optimizations
> - KV pages of 16, 32, or 64 tokens, matching paged-attention serving systems. Each page stores a compact header containing scale exponents, residual radius, stage offsets, and optional sensitivity metadata
> - Stage-0 prefix $\in$ fast memory, directly consumed by the attention kernel. 
> 	- For GPU: start with FP8 or INT8, then add a 4-bit/block-scaled path for newer targets
> 	- For `Serif`/FPGA:  packed INT4 or INT2 and consumed by fused DP4A-style dot products
> - Each refinement stage should encode the residual of the previous stage: $$r^{(m)} = x - \hat{x}^{(m-1)}, $$ then  $$\hat{x}^{(m)} = \hat{x}^{(m-1)} + \Delta^{(m)},$$ which makes every prefix useful and lets stage $m+1$ strictly tighten radius $\rho_m$ 
> - $\forall p, h, g, m$ (page, head, group, and stage) store or derive $\rho_{p,h,g}^{(m)}$
> 	- Deterministic version uses quantization-bin geometry
> 	- Empirical high-confidence version uses calibration quantities
> 		- Report both
> - Refine keys first and only values whose certified softmax mass is non-negligible to reduce refinement traffic
> - For each page, compute a score upper bound and only pages with competitive upper bounds enter the candidate set
> 	- Safer than ordinary top-$k$ sparse attention because a page can be dropped only when its upper bound proves it cannot carry enough mass
> - Implement all interval calculations with log-sum-exp stabilization:
> 	- For ignored pages $I$ and active pages $A$, $$\log N_I=\log\sum_{i\in I}e^{U_i},\qquad \log D_a=\log\sum_{i\in A}e^{L_j},$$ then the omitted-mass bound is $$\delta_I=\frac{e^{\log N_I}}{e^{\log D_A}+e^{\log N_I}}$$ which avoids numerical instability @ long context
> - Keep the recent window and known attention-sink tokens at a higher stage by default
> - If a page was refined at token $t$, keep that refinement resident for a short horizon $H$ 
> 	- Locality = adjacent decode steps often inspect overlapping pages
> - Use refine when the error bound exceeds $\epsilon_{\text{high}}$ threshold and evict when it falls below $\epsilon_{\text{low}}$ threshold
> - Multiple requests can need refinement from the same layer/head/page stage in serving so group fetches into one memory transaction (stage 1 + 2 $\in$ CXL memory?)
> - Within each KV page, quantize token streams with causal residual feedback
> 	- Shapes value error into high token-frequency modes, which smooth attention distributions suppress
> 	- For keys, run the shaping in pre-RoPE or RoPE-demodulated coordinates to avoid treating deterministic rotary phase as signal variation
> - Use no shaping for spike/retrieval heads, first-order shaping for smooth heads, and possibly second-order shaping for very smooth summarization heads
> 	- Instead choice should come from a measured token-axis spectral statistic
> - Allocate higher stage-0 precision to heads with high certified sensitivity $S_h \approx |A_h D_{p_h}|_F^2,$
> 	- $A_h$ is the empirical read operator and **$D_{p_h}$ is the shaping filter from our write-time decoder**
> 	- Low-sensitivity heads get fewer resident bits
> - For MoE models, store coarse gate projections or gate-weight prefixes in HBM then compute intervals for gate logits
> 	- If top-$k$ experts are certified, do not fetch deeper gate precision
> 	- If the selected experts are also stored progressively, fetch only their needed weight refinements
> - After computing a low-stage forward pass, bound the final hidden-state error
> 	- If the top vocabulary logit remains separated under that bound, accept the token without further refinement
> 	- If not, refine (second stopping rule)
> - Give two modes
> 	- **Exact**, use worst-case radii
> 	- **Probabilistic**, use calibrated high-confidence radii, e.g., 99.9% per-page or per-head bounds
> - Implement a custom attention module that plugs into HuggingFace/vLLM-style decoding
> 	- Start with a PyTorch/Triton prototype, then move the bottleneck pieces into custom CUDA/Triton kernels
> - The prefix path should produce the provisional attention result quickly
> 	- Refinement correction should run only on selected pages and can be accumulated as a correction term
> - Every experiment should **report HBM bytes/token, host bytes/token, refinement bytes/token, and useful arithmetic/byte**

---

## Core Theorems
1. [[#^de51de|Progressive code certificate]]
2. [[#^b8a0f7|Certified softmax-tail bound]]
3. [[#^287a53|Certified attention-output bound]]
4. [[#^94c4cd|MoE exact top-k routing certificate]]
5. [[#^20c7b9|Expected refinement theorem]]
6. [[#^6f9ac8|Expected byte theorem]]
7. [[#^f72505|Read-aware bit allocation]]
8. [[#^5c64b2|Write-time encoder value-path theorem]]
9. [[#^016da2|Write-time encoder key-path first-order theorem]]
10. [[#^fa9848|Read-aware rate-distortion lower bound]]
11. [[#^a10c24|Output-token certification theorem]]
---
> [!question] Progressive code certificate
> For each block $x_i$, stage $m$ reconstruction $\hat{x}_i^{(m)}$ satisfies $$||x_i-\hat{x}_i^{(m)}||_2 \le \rho_i^{(m)}$$ and for attention keys $$s_i = \frac{q^\top k_i}{\sqrt d},\qquad \hat{s}_i^{(m)} = \frac{q^\top \hat{k}_i^{(m)}}{\sqrt d},$$ then $$|s_i-\hat{s}_i^{(m)}|\le\frac{||q||_2}{\sqrt d}\rho_{K,i}^{(m)}=\alpha_i^{(m)},$$ so every score has a certified interval $$L_i^{(m)}=\hat{s}_i^{(m)}-\alpha_i^{(m)},\qquad U_i^{(m)}=\hat{s}_i^{(m)}+\alpha_i^{(m)},$$ our main primitive.

^de51de

---
> [!question] Certified softmax-tail bound
> Let $A$ be the set of active/refined tokens or pages and $I$ the ignored/unrefined set. For true softmax probabilities $$p_i=\frac{e^{s_i}}{\sum_j e^{s_j}},$$ the total true mass on ignored tokens satisfies $$\sum_{i\in I}p_i\le\delta_I=\frac{\sum_{i\in I}e^{U_i}}{\sum_{j\in A}e^{L_j}+\sum_{i\in I}e^{U_j}},$$ which is exact and deterministic given the intervals.

^b8a0f7

---
> [!question] Certified attention-output bound
> Assume $|v_i|_2\le V$. If the ignored set has mass at most $\delta_I$, then replacing the full attention output with attention over only $A$ incurs $$||o-o_A||_2\le 2V\delta_I.$$ If active-set scores are also approximate with $$\epsilon_A=\max_{i\in A}|s_i-\hat{s}_i|,$$ then softmax perturbation contributes at most $2V\tanh(\epsilon_A).$ If values are reconstructed with $$||v_i-\hat v_i||_2 \le \eta_i,$$ then value error contributes $\sum_{i\in A}\hat p_i\eta_i.$  A clean full certificate is therefore $$||o-\hat o||_2\le 2V\delta_I + 2V\tanh(\epsilon_A)+\sum_{i\in A}\hat p_i\eta_i,$$ which gives the runtime a direct stopping condition $$2V\delta_I + 2V\tanh(\epsilon_A)+\sum_{i\in A}\hat p_i\eta_i\le\varepsilon_{\text{attn}}.$$

^287a53

---
> [!question] MoE exact top-$k$ routing certificate
> For expert gate logits $g_e$ suppose PRISM computes intervals $g_e\in[L_e,U_e]$ and let $S$ be the set of $k$ experts with the largest lower bounds. If $\min_{e\in S} L_e = \max_{f\notin S}U_f,$ then $S$ is exactly the true top-$k$ expert set, and no further gate refinement is needed.

^94c4cd

---
> [!question] Expected refinement theorem
> Let $\tau$ be a live decision threshold: an attention retention threshold, a page-selection boundary, or an MoE top-$k$ boundary, then define the true margin $$\Delta_i = |s_i-\tau|.$$ Then at stage $m$, block $i$ is uncertain only if its interval overlaps the threshold and since $|s_i-\hat{s}_i^{(m)}|\le \alpha_i^{(m)},$ uncertainty implies approximately $\Delta_i \le 2\alpha_i^{(m)}.$ If the margin density near zero is bounded by $G$, then $$\Pr(\text{block } i \text{ remains uncertain at stage } m)\le 4G\alpha_i^{(m)}.$$ Using $$\alpha_i^{(m)}=\frac{||q||_2}{\sqrt d}\rho_i^{(m)},$$ we get $$\Pr(\text{uncertain})\le\frac{4G||q||_2}{sqrt d}\rho_i^{(m)},$$ and if the progressive codec gives $\rho_i^{(m)}\le C_i 2^{-R_m}$ then $$\Pr(\text{uncertain})\le\frac{4G||q||_2C_i}{\sqrt d}2^{-R_m},$$ meaning expected refinement traffic decays exponentially with prefix rate provided the model has non-pathological decision margins.

^20c7b9

---
> [!question] Expected byte theorem
> Let $c_m$ be the byte cost of fetching refinement stage $m$ and let $\pi_m$ be the probability that a block survives to stage $m$. Then expected bytes/token are $$\mathbb{E}[\text{bytes/token}]=\sum_i c_{i,0}+\sum_i\sum_{m\ge 1}c_{i,m}\pi_{i,m},$$ and then from the expected refinement bound, $$\pi_{i,m}\lesssim\frac{4G_i|q|_2C_i}{\sqrt d}2^{-R,m}.$$ Therefore, PRISM's deeper-memory is controlled by margin density and not solely context length.

^6f9ac8

---
> [!question] Read-aware bit allocation
> Suppose head/page/group $h$ has sensitivity $S_h$, dynamic range $R_h$ and bit allocation $b_h$. High-rate quantization gives approx. $$D_h(b_h)\approx cR_h^2s_H2^{-2b_h}$$ with a global bit budget $\sum_h n_h b_h\le B,$ so the Lagrangian optimum is $$b_h^\star=\frac12\log_2\left(\frac{cR_h^2S_h}{\lambda n_h}\right)_+.$$ More bits go to high-sensitivity, high-dynamic-range groups.

^f72505

---
> [!question] Write-time encoder value-path theorem
> For a scalar token stream $x_t$, define causal residual-feedback quantization $$u_t=x_t+r_{t-1},$$ $$\hat{x}_t=Q(u_t),$$ $$r_t=u_t-\hat{x}_t.$$ The quantization error is $$\varepsilon_t=\hat{x}_t-x_t=r_{t-1}-r_t.$$ $\forall$ future attention read $y=\sum_{t=1}^T a_tx_t,$ the readout error is $$\hat y-y=\sum_{t=1}^T a_t(r_{t-1}-r_t),$$ and since $r_0=0$ this telescopes to $$\hat y-y=\sum_{t=1}^{T-1}(a_{t+1}-a_t)r_t-a_Tr_T.$$ If $|r_t|\le\Delta/2$, then $$|\hat y-y|\le\frac{\Delta}{2}\left(\operatorname{TV}(a)+|a_T|\right),$$ where $\operatorname{TV}(a)=\sum_{t=1}^{T-1}|a_{t+1}-a_t|.$ 
> For values, $a_t=p_{\tau,t}$ so $$|\delta o_{\tau,c}^{(V)}|\le\frac{\Delta_c}{2}\left(\operatorname{TV}(p_\tau)+p_{\tau,T}\right),$$ proving the write-time encoder error is small when attention is smooth over token order.

^5c64b2

---

> [!question] Write-time encoder key-path first-order theorem
> For keys, $$s_{\tau,t}=\frac{q_\tau^\top k_t}{\sqrt d},\qquad p_\tau=\operatorname{softmax}(s_\tau),\qquad o_\tau=\sum_t p_{\tau,t}v_t.$$ A first-order perturbation gives $$\delta o_\tau \approx \sum_t p_{\tau,t}\delta v_t + \sum_t p_{\tau,t}(v_t-o_\tau)\frac{q_\tau^\top\delta k_t}{\sqrt d}.$$ The value path is weighted by $p_{\tau,t}$ and the key path by $p_{\tau,t}(v_t-o_\tau)q_{\tau,j}$; applying the same telescoping argument coordinatewise shows that shaped key error is controlled by the token-axis variation of $$b_{\tau,t}=p_{\tau,t}(v_t-o_\tau).$$

^016da2

---
> [!question] Read-aware rate-distortion lower bound
> For a Gaussian block $x\sim \mathcal N(0,\Sigma),$ and a downstream read operator $A$, the distortion is $$d_A(x,\hat{x})=|A(x-\hat{x})|_2^2.$$ Let $\lambda_i$ be the eigenvalues of $$\Sigma^{1/2}A^\top A\Sigma^{1/2}.$$ Then the Gaussian weighted rate-distortion function is $$R(D)=\frac12\sum_i\log_+\left(\frac{\lambda_i}{\theta}\right)$$ with $$D=\sum_i\min(\lambda_i,\theta).$$ Optimal quantization should spend bits on modes that the future read operator actually observes.

^fa9848

---
> [!question] Output-token certification theorem
> Let the approximate final hidden state be $\hat h$, with $$|h-\hat h|_2\le \epsilon_h.$$ For vocabulary row $w_y$, the logit error satisfies $$|w_y^\top h-w_y^\top \hat h|\le||w_y||_2\epsilon_h.$$ Define logit intervals $z_y\in[\hat z_y-||w_y||_2\epsilon_h, \hat z_y+||w_y||_2\epsilon_h].$ If token $a$ satisfies $$\hat z_a-||w_a||_2\epsilon_h>\max_{y\ne a}\left(\hat z_y+||w_y||_2\epsilon_h\right),$$ then the next token is exactly the same as the full-precision model’s greedy output.

^a10c24

---
> [!warning] Literature distinction
>`TurboQuant` (and really `RaBitQ`) gives near-optimal static vector quantization through random rotation, scalar optimal quantizers, and a QJL residual for unbiased inner products. `RaBitQ` gives an unbiased distance estimator and sharp $O(1/\sqrt D)$ error bound for ANN-style high-dimensional vector search. PRISM optimizes certified runtime decisions instead of vector geometry.

---
## Experimental plan

> [!success] Margin audit
>  Determining experimental quality empirically. For long-context dense models, measured attention page margins, omitted softmax mass, and how many pages remain ambiguous after 2-bit, 3-bit, and 4-bit prefixes. For MoE models, measured top-$k$ gate margins and how often the selected experts are certified from coarse prefixes. Most pages/experts were far from the decision boundary, so PRISM is a viable technique.

- [ ] **Validate certificates**. For attention, plot true omitted mass against the certified $\delta_I$. For output error, plot true attention-output error against the theoretical bound. For MoE, plot exact top-$k$ agreement under interval certification. ==(*Weeks 1, 2*)==
- [ ] **GPU prototype.** Implement PRISM-KV in PyTorch/Triton/CUDA with stage-0 resident prefixes and stage-1 residual fetch. Report tokens/sec, inter-token latency, HBM bytes/token, host bytes/token, and accuracy. Compare against full FP16/BF16 KV, static INT8/FP8 KV, static INT4 with dequantization, KIVI-style KV quantization, TurboQuant-like static quantization, and dequant-wall measurements from `Serif` ==(*Week 2*)==
- [ ] **Heterogenous memory.** Keep stage 0 in GPU HBM, stage 1 in pinned host memory, and optionally stage 2 in CPU memory or an emulated CXL tier. Ideally would get the same or near-same quality while fetching only a small fraction of deeper bytes per token. ==(*Week 2,3*)==
- [ ] **MoE.** Using model where expert memory pressure matters and comparing normal expert loading/routing against PRISM-Gate. Report the fraction of tokens whose top-$k$ experts are certified from the coarse prefix, the number of expert refinements per token, and end-to-end throughput. ==*(Week 3)*==
- [ ] **Failure analysis.** Construct hard prompts with needle retrieval, sharp attention spikes, repeated keys, adversarially similar MoE gate logits, and low-margin output distributions. PRISM should degrade gracefully such that when the certificate cannot prove safety, it fetches more bits and falls back toward full precision. ==*(Week 3, 4)*==