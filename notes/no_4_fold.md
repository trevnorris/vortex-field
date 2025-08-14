# executive summary (what changes & why)

**Root cause of the “×4”**

* a normalization mistake in the 4D→3D **kernel prefactor** (used $1/(2\pi)$ where $1/(4\pi)$ belongs),
* **double-counting** the $w=0$ neighborhood as a “direct” term in addition to integrating over $w$, and
* attributing a nonzero loop contribution to a **continuity/“drainage” potential** (which has zero circulation).

**Correct statement (projection invariance)**
For a codimension-2 sheet $\Sigma\subset\mathbb{R}^4$ with sheet strength $\Gamma$, and any small loop $\gamma\subset\Pi=\{w=0\}$ linking $\Sigma\cap\Pi$ once,

$$
\oint_{\gamma}\mathbf v(\cdot,0)\cdot d\boldsymbol\ell
\;=\;
\Gamma\;+\;O\big((\xi/\rho)^2+(\kappa\rho)^2\big).
$$

In a half-space split, **each half-space contributes $\Gamma/2$ to the *circulation***; there is **no separate “direct”** or “drainage” add-on.

**GR weak-field/GEM**
The coefficient $16\pi G/c^2$ in the vector Poisson equation comes from **linearized GR with trace reversal** and the identification of $\mathbf A_g$ with $\bar h_{0i}$ in Lorenz gauge. It **does not require** a geometric $4\times$ factor. Keep that coefficient; remove any “4× decomposition” narrative.

---

# global “find & remove/replace” checklist

Search the .tex for the following (case-insensitive):

* `4-fold`, `fourfold`, `4\\times`, `four-fold`, `four fold`

  * **Action:** remove claims of “exact fourfold enhancement”; replace with “projection invariance of circulation”.
* “direct + upper + lower + drainage = 4Γ”

  * **Action:** remove “direct” and “drainage” as independent adders; keep only “upper” and “lower” half-spaces, each $\Gamma/2$, summing to $\Gamma$.
* any text asserting “continuity/drainage induces a swirl that adds Γ”

  * **Action:** replace with: “the continuity-driven potential part does **not** contribute to loop circulation.”
* any place that algebraically factors $16\pi G/c^2$ as $4\times4\times(\pi G/c^2)$

  * **Action:** delete that factoring and the surrounding justification; insert a short linearized-GR derivation (provided below).

---

# exact latex replacements & inserts

## 0) title & abstract

**Title (if it mentions fourfold):**
Replace with something like:

```
\title{Projection Invariance of Circulation from 4D$\to$3D Vortex Slicing}
```

**Abstract (surgical edits):**

* Replace “…we prove an exact \emph{four-fold} enhancement…” →
  “…we prove an exact \emph{projection invariance}: the measured 3D circulation equals the intrinsic sheet strength \$\Gamma\$ under clear geometric/regularity assumptions.”
* Replace the “four equal contributions” sentence with:
  “A half-space split shows each side (\$w>0\$ and \$w<0\$) contributes \$\Gamma/2\$ to the loop circulation on the slice; potential (‘drainage’) parts do not affect loop circulation.”

## 1) postulates section

**Postulate 5 (currently: “Quantized vortices with 4-fold projection”)**
Replace the header & body with:

```
\textbf{Postulate 5 (Projection invariance of circulation).}
For a codimension-2 vortex sheet $\Sigma\subset\mathbb{R}^4$ with sheet strength $\Gamma$, the loop circulation measured on the slice $\Pi=\{w=0\}$ by any small loop $\gamma$ linking $\Sigma\cap\Pi$ once equals $\Gamma$ in the thin/flat limit ($\xi/\rho\to 0$), with corrections $O((\xi/\rho)^2+(\kappa\rho)^2)$.
```

**Postulate 6 (Discrete vortex projection)**
Leave it, but **remove** any language that says the projection multiplies circulation; it may still discuss discrete intersection/occupancy if you like (as *sources*, not as a geometric factor).

## 2) main theorem & lemmas

**Delete** the “Four-fold projection” theorem and its four lemmas.
**Insert** the following compact theorem + two lemmas:

```
\begin{theorem}[Projection invariance]
Let $\Sigma\subset\mathbb{R}^4$ be a smooth codimension-2 vortex sheet of strength $\Gamma$, $\Pi=\{w=0\}$, and $\gamma_\rho\subset\Pi$ a small loop linking $\mathcal{C}_0=\Sigma\cap\Pi$ once. In the thin/flat limit $\xi/\rho\to 0$ with curvature radius $\kappa^{-1}\gg \rho$, 
\[
\oint_{\gamma_\rho}\mathbf v(\cdot,0)\cdot d\boldsymbol\ell
= \Gamma + O\!\big((\xi/\rho)^2+(\kappa\rho)^2\big).
\]
\end{theorem}

\begin{lemma}[Half-space split]
With the correct 4D $\to$ 3D kernel normalization,
\[
v_\theta(\rho)
=\frac{\Gamma}{4\pi\rho}\!\int_{-\infty}^{\infty}\!\frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}
=\frac{\Gamma}{2\pi\rho}.
\]
Consequently, each half-space contributes
\[
\oint_{\gamma_\rho}\mathbf v^{(\pm)}\cdot d\boldsymbol\ell = \frac{\Gamma}{2}.
\]
\end{lemma}

\begin{lemma}[No potential contribution]
On the slice $\Pi$, Helmholtz-decompose $\mathbf v=\nabla\phi + \nabla\times \mathbf A$. The continuity-driven “drainage’’ part is potential ($\nabla\phi$), hence
\(
\oint_{\gamma_\rho}\nabla\phi\cdot d\boldsymbol\ell = 0.
\)
\end{lemma}
```

## 3) kernel section (fix normalization)

Replace the kernel/prefactor block with:

```
\paragraph{Projected kernel with correct normalization.}
In $\mathbb{R}^4$ the Green's function is $G_4(x)=-1/(4\pi^2|x|^2)$, so velocity contributions decay like $|x|^{-3}$. For the azimuthal component on the slice $\Pi$ induced by material at offset $w$,
\[
K(\rho,w)=\frac{\rho^2}{(\rho^2+w^2)^{3/2}},\qquad
v_\theta(\rho)
=\frac{\Gamma}{4\pi\rho}\!\int_{-\infty}^{\infty}\!K(\rho,w)\,dw.
\]
The identities
\(
\int_0^\infty \frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}=1
\)
and
\(
\int_{-\infty}^{\infty} \frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}=2
\)
give $v_\theta(\rho)=\Gamma/(2\pi\rho)$ and each half-space yields $\Gamma/2$ circulation.
```

Remove any statement that integrates $1/(\rho^2+w^2)$ or uses a $1/(2\pi)$ prefactor.

## 4) drainage/continuity paragraph

Replace the “drainage = Γ” claim with:

```
\paragraph{Continuity and potential flow on $\Pi$.}
Axial flux through a pillbox straddling $\Pi$ induces a potential adjustment $\nabla\phi$ on the slice to satisfy mass conservation. Being gradient, this part does not contribute to loop circulation: $\oint_\gamma\nabla\phi\cdot d\boldsymbol\ell=0$.
```

## 5) error bounds (tighten justification)

Replace the “two bounded moments” line with the convolution/Taylor estimate:

```
\paragraph{Finite thickness (mollified profile).}
Let $K_\rho(w)=\rho^2(\rho^2+w^2)^{-3/2}$ and let $\eta_\xi$ be an even mollifier with unit mass and second moment $\mu_2=O(\xi^2)$. Then
\[
\Big|\!\int (\eta_\xi*K_\rho)(w)\,dw - \int K_\rho(w)\,dw\Big|
\le \tfrac{\mu_2}{2}\,\|\partial_w^2 K_\rho\|_{L^1}
= O\big((\xi/\rho)^2\big).
\]
```

Keep the curvature $O((\kappa\rho)^2)$ statement as-is.

## 6) “relation to prior work” sentence

Change “this paper fills the gap by proving an exact fourfold law” to “this paper gives a compact proof of **projection invariance** (circulation equals sheet strength) in 4D→3D slicing with explicit kernel normalization and error control.”

## 7) vector (GEM) sector in GR (where you previously used 4×)

Where you present the vector equation with $16\pi G/c^2$, **delete** any factoring/explanation that allocates a “4” to projection. Insert this short derivation box:

```
\paragraph{Linearized GR normalization of the vector sector.}
In Lorenz gauge, trace-reversed metric $\bar h_{\mu\nu}=h_{\mu\nu}-\tfrac12\eta_{\mu\nu}h$ obeys
\[
\square\,\bar h_{\mu\nu}=-\frac{16\pi G}{c^4}\,T_{\mu\nu}.
\]
With slow motion $T_{00}\approx \rho c^2$, $T_{0i}\approx \rho c\,v_i=J_i c$, define
\[
\Phi_g=\frac{c^2}{2}\,\bar h_{00},\qquad
\mathbf A_g=-\,\frac{c^2}{4}\,\bar{\mathbf h}_{0i}\,\hat{\mathbf e}_i
\]
so that the Lorenz condition becomes the standard GEM gauge. Then
\[
\nabla^2\Phi_g=4\pi G\,\rho,\qquad
\nabla^2\mathbf A_g=-\,\frac{16\pi G}{c^2}\,\mathbf J.
\]
This fixes the $16\pi G/c^2$ coefficient independently of any projection factor.
```

If you prefer to keep a one-knob calibration, you may add a dimensionless $\kappa$ as a single multiplier and set $\kappa=1$ by Lense–Thirring phenomenology; be explicit that $\kappa$ is fixed by GR, not by projection.

---

# things to delete outright

* The “Four equal contributions ($\Gamma+\Gamma+\Gamma+\Gamma$)” figure/text and any equation summing to $4\Gamma$.
* The lemma claiming “drainage” contributes $\Gamma$ to loop circulation.
* Any derivation that uses $\int_0^\infty dw/(\rho^2+w^2)=\pi/(2\rho)$ for the projected kernel. (Wrong kernel/power for the 4D slice problem.)
* Any factorization line explaining $16\pi G/c^2$ as $4\times4\times(\pi G/c^2)$.

---

# what to keep (unchanged)

* Dimensional analysis/units checks (they remain correct once the prefactor is fixed).
* Thin/flat limit assumptions; curvature and thickness corrections $O((\kappa\rho)^2)+O((\xi/\rho)^2)$.
* Worked geometries (“straight pierce”, “tilted disk”) — just ensure they now illustrate **$\Gamma/2+\Gamma/2=\Gamma$**, not $4\Gamma$.

---

# minimal numeric verification plan (to include in appendix)

**Goal:** verify half-spaces $\rightarrow \Gamma/2$ each; total $\rightarrow \Gamma$.

1. Choose an even mollifier $\sigma_\xi(w)$ with unit mass (e.g., Gaussian).
2. Compute

$$
I_\pm(\rho,\xi)=\int_{0}^{\infty}\sigma_\xi(w)\,\frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}
\quad\text{and}\quad
I_-(\rho,\xi)=I_+(\rho,\xi).
$$

3. Set $v_\theta^{(\pm)}(\rho)=\frac{\Gamma}{4\pi\rho}\,I_\pm$, then $\mathcal C_\pm=2\pi\rho\,v_\theta^{(\pm)}$.
4. Check $\mathcal C_+\approx \mathcal C_-\approx \Gamma/2$ and $\mathcal C_++\mathcal C_- \approx \Gamma$ as $\xi/\rho\to 0$.
5. Verify that adding any gradient field on $\Pi$ does not change $\oint v\cdot d\ell$.

(If you want, I can provide a short python/sympy block you can drop into an appendix.)

---

# diff-style snippets you can paste

### theorem block (drop-in)

```
\begin{theorem}[Projection invariance]
Under the thin/flat assumptions,
\[
\oint_{\gamma_\rho}\mathbf v(\cdot,0)\cdot d\boldsymbol\ell
= \Gamma + O\!\big((\xi/\rho)^2+(\kappa\rho)^2\big).
\]
\end{theorem}
```

### kernel identity (drop-in)

```
\[
\int_0^\infty \frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}=1,
\qquad
\int_{-\infty}^{\infty} \frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}=2.
\]
\[
v_\theta(\rho)
=\frac{\Gamma}{4\pi\rho}\!\int_{-\infty}^{\infty}\frac{\rho^2\,dw}{(\rho^2+w^2)^{3/2}}
=\frac{\Gamma}{2\pi\rho}.
\]
```

### forms/Stokes framing (clean proof sketch)

```
Let $u$ be the velocity 1-form on $\mathbb{R}^4$, $\Omega=du$ the vorticity 2-form (a de Rham current supported on $\Sigma$). Let $S\subset\Pi$ be a small ribbon with $\partial S=\gamma_\rho$ and let $T$ be a small cap linking $\Sigma$ once. Since $du=\Omega$ and $d\Omega=0$ away from $\Sigma$,
\[
\oint_{\gamma_\rho}u = \int_S du = -\int_T du = -\int_T \Omega = \Gamma,
\]
by the definition of sheet strength via transverse flux. Curvature/thickness give the stated $O((\kappa\rho)^2)+O((\xi/\rho)^2)$ corrections.
```

### GEM normalization insert

(use the derivation box in §7 above)

---

# risks to scan for after edits

* Any downstream formula that literally uses $4\Gamma$ as an input; replace with $\Gamma$.
* Any discussion of experiment/simulation “validating” a 4×. Recast as validating the **invariance** (sum = $\Gamma$) or remove.
* Figures/tables labeling “×4” — update captions.

---

# one-paragraph “why our opinion changed” (optional preface for the paper)

> A previous draft used an exact “four-fold enhancement” derived from a projected kernel formula. Revisiting the 4D Green’s normalization and avoiding double counting at the $w=0$ slice shows that circulation on the 3D slice equals the sheet strength $\Gamma$ (projection invariance). A half-space split contributes $\Gamma/2$ per side; potential (continuity-driven) terms produce no loop circulation. The weak-field GR vector coefficient $16\pi G/c^2$ is then obtained from linearized GR via the standard trace-reversed formulation, independently of any projection factor.

---

# deliverables you’ll hand me in the next session

* The .tex with:

  * all “4-fold” occurrences highlighted or already removed,
  * the theorem/lemmas replaced by the projection-invariance version,
  * the kernel block corrected,
  * the GEM normalization paragraph inserted,
  * any figure/table captions updated.
* If you’d like, a short appendix code cell (python/sympy or numpy) to numerically verify the half-space split; I can write/clean that on the spot.

that’s the whole repair kit. once you paste this into the doc, the narrative is actually cleaner and more defensible, and your GR mapping no longer relies on a fragile artifact.

