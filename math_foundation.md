# Math Foundations: Eigenvectors, Orthogonality, and Fourier Projection

A step-by-step build-up of the linear algebra ideas that underlie Fourier analysis and signal processing.

---

## Section 1 — Eigenvectors and Eigenvalues

### Step 1: A vector is just an arrow

A vector in 2D is a pair of numbers like `v = [2, 1]`. Picture it as an arrow from the origin to the point (2, 1). It has a **direction** and a **length**.

### Step 2: A matrix is a "transformation machine"

A matrix takes a vector and turns it into another vector. For example:

```
A = [2  0]
    [0  3]
```

If we feed in `v = [1, 1]`:

```
A · v = [2·1 + 0·1]   = [2]
        [0·1 + 3·1]     [3]
```

So the matrix **stretched** the arrow — it changed both its direction and length.

In general, when a matrix acts on a vector, it can:
- stretch it
- shrink it
- rotate it
- flip it
- or any combination

### Step 3: The special vectors — eigenvectors

Here's the magic question:
> Are there any special vectors that, when the matrix acts on them, **don't change direction** — only get stretched or shrunk?

If yes, those are called **eigenvectors**. ("Eigen" is German for "own" or "characteristic".)

Formally:

```
A · v = λ · v
```

In words: applying matrix A to vector v gives you the SAME vector v back, just scaled by some number λ (lambda). That number λ is the **eigenvalue**.

- `v` = the eigenvector (direction that survives unchanged)
- `λ` = the eigenvalue (how much it got stretched, shrunk, or flipped)

### Step 4: Concrete example

Take the matrix from before:
```
A = [2  0]
    [0  3]
```

Try `v = [1, 0]`:
```
A · v = [2·1 + 0·0]   = [2]   = 2 · [1]
        [0·1 + 3·0]     [0]       [0]
```
Same direction, scaled by 2. So `[1, 0]` is an eigenvector, eigenvalue = 2.

Try `v = [0, 1]`:
```
A · v = [0]   = 3 · [0]
        [3]       [1]
```
Same direction, scaled by 3. So `[0, 1]` is an eigenvector, eigenvalue = 3.

### Step 5: Why we care

Eigenvectors reveal the **natural axes** of a transformation — the directions in which a complicated transformation acts like a simple stretch. If you express problems in terms of eigenvectors, complicated matrix operations collapse into simple multiplications. This is the engine behind PCA, vibration analysis, quantum mechanics, Google PageRank, and — as we'll see — Fourier analysis.

---

## Section 2 — Vector Orthogonality

### Step 1: The dot product

Given two vectors `a = [a₁, a₂]` and `b = [b₁, b₂]`, their **dot product** is:

```
a · b = a₁·b₁ + a₂·b₂
```

Geometrically:
```
a · b = |a| · |b| · cos(θ)
```
where θ is the angle between them.

### Step 2: Orthogonal = perpendicular

If the angle θ = 90°, then cos(90°) = 0, so:

```
a · b = 0   ⟺   a ⊥ b
```

Two vectors are **orthogonal** when their dot product is zero. Example: `[1, 0]` and `[0, 1]` — clearly perpendicular, and `1·0 + 0·1 = 0`. ✓

### Step 3: Why orthogonality is powerful

If your basis vectors are orthogonal, you can decompose any vector cleanly, without the components "interfering":

```
v = c₁·e₁ + c₂·e₂ + c₃·e₃ + ...
```

Each coefficient `cᵢ` is found just by taking the dot product `v · eᵢ`. No solving systems of equations, no leakage between components. Each axis is **independent**.

Think of mixing paint: if your basis colors are truly independent (orthogonal), you can extract the exact amount of red, blue, green from any mixture by a single operation per color. If they weren't independent, extracting "red" would also pull some blue with it.

#### Why does `cᵢ = v · eᵢ` actually work? — the derivation

This is the key trick, so let's see it carefully.

**Step A — set up.** We know `v` can be written as some combination of the basis vectors with unknown coefficients:

```
v = c₁·e₁ + c₂·e₂ + c₃·e₃ + ...
```

The question: how do we find `c₁`?

**Step B — dot both sides with `e₁`.**

```
v · e₁ = (c₁·e₁ + c₂·e₂ + c₃·e₃ + ...) · e₁
```

Distribute the dot product:

```
v · e₁ = c₁·(e₁·e₁) + c₂·(e₂·e₁) + c₃·(e₃·e₁) + ...
```

**Step C — use orthogonality to kill every term except one.**

By definition, an orthogonal basis satisfies `eᵢ · eⱼ = 0` whenever `i ≠ j`. So:

```
v · e₁ = c₁·(e₁·e₁) + c₂·  0   + c₃·  0   + ...
                            ↑           ↑
                  killed because    killed because
                  e₂ ⊥ e₁          e₃ ⊥ e₁
```

Only the first term survives:

```
v · e₁ = c₁ · (e₁ · e₁)
```

**Step D — if `e₁` has unit length, `e₁·e₁ = 1`, so:**

```
c₁ = v · e₁
```

Done. And the same trick gives `c₂ = v · e₂`, `c₃ = v · e₃`, etc. — each coefficient is one dot product.

#### Concrete example — orthogonal basis (it works)

Let `v = [3, 4]`, basis `e₁ = [1, 0]`, `e₂ = [0, 1]`.

```
c₁ = v · e₁ = 3·1 + 4·0 = 3   ✓
c₂ = v · e₂ = 3·0 + 4·1 = 4   ✓
```

So `v = 3·e₁ + 4·e₂`. Sanity check: `3·[1,0] + 4·[0,1] = [3,4]`. ✓

#### Counter-example — NON-orthogonal basis (it breaks)

Take a basis that is **not** orthogonal:

```
f₁ = [1, 0]
f₂ = [1, 1]
```

Check: `f₁ · f₂ = 1·1 + 0·1 = 1` ≠ 0 → not orthogonal.

The **true** decomposition of `v = [3, 4]` into this basis comes from solving the system:

```
c₁·[1,0] + c₂·[1,1] = [3,4]
  → c₁ + c₂ = 3
  →       c₂ = 4
  → c₂ = 4, c₁ = −1
```

So the correct coefficients are `c₁ = −1`, `c₂ = 4`.

But what does the dot-product shortcut give?

```
v · f₁ = 3·1 + 4·0 = 3
```

It says `c₁ = 3`. **But the real answer is c₁ = −1.** The shortcut lied.

Why? Redo Step C without orthogonality:

```
v · f₁ = c₁·(f₁·f₁) + c₂·(f₂·f₁)
       = c₁·1       + c₂·1
       = c₁ + c₂
```

The `c₂` term **didn't die** — because `f₂` wasn't perpendicular to `f₁`. The dot product mixes contributions from `c₁` AND `c₂` together. You can't read off `c₁` cleanly — you'd have to solve a system of equations to untangle them.

#### The big picture

| Basis is...     | Effect of `v · eᵢ`                              |
|-----------------|-------------------------------------------------|
| Orthogonal      | Only the `cᵢ` term survives → clean readout     |
| Not orthogonal  | All terms leak in → tangled mess                |

Orthogonality is exactly the property that makes coefficients readable by a single, independent measurement — one dot product per direction, no cross-talk. That's why we'll need sin/cos at different frequencies to be mutually orthogonal in Section 4: it's what lets each frequency's coefficient be one clean integral, with no leakage from other frequencies.

---

## Section 3 — Fourier Series, the Fourier Transform, and Euler's Formula

Before we use orthogonality on sines and cosines, let's see the bigger picture: what is a Fourier series, what is the Fourier transform, and why do we suddenly need `e^(iωt)` instead of plain sines and cosines?

### Step 1: The time domain

When you record a signal — a voice, a heartbeat, a stock price — you collect a value at each moment in time. Plotted, the x-axis is time `t` and the y-axis is the signal value `f(t)`. This is the **time domain**.

```
f(t) ▲
     │      ╭─╮      ╭─╮
     │     ╱   ╲    ╱   ╲
     │────╯     ╲──╯     ╲────→  t
     │           ╲╱
```

It's how we naturally observe signals.

### Step 2: The frequency domain — a different lens on the same signal

Same signal, different question:
> "What pure tones (sinusoids) is this signal made of, and how loud is each?"

Answer: a list of amplitudes, one per frequency.

```
|F(ω)| ▲
       │           ●  ← lots of frequency-3 content
       │
       │     ●        ← some frequency-1 content
       │
       │  ●        ●     ← bits of other frequencies
       └──────────────────→ ω
         0  1  2  3  4  5
```

This list IS the **frequency domain** representation of the same signal.

Two views, one signal:

| Time domain | Frequency domain |
|-------------|------------------|
| `f(t)` — how the signal changes over time | `F(ω)` — how much of each frequency the signal contains |
| Natural for recording / playback | Natural for filtering, compression, tone identification |
| "What is happening at this moment?" | "What is the signal made of?" |

Analogy: time domain = a street map of a city; frequency domain = an ingredient list (steel, glass, concrete). Same city, different language. The **Fourier transform** is the translator between the two.

### Step 3: Fourier series — periodic signals as sums of sines and cosines

Joseph Fourier (1822): **any reasonable periodic signal** can be written as a sum of sines and cosines:

```
              ∞
f(t) = a₀/2 + Σ [aₙ·cos(nω₀t) + bₙ·sin(nω₀t)]
             n=1

where ω₀ = 2π/T  (T is the period)
```

The pieces:
- `a₀/2` — the constant offset (DC level)
- `n = 1`: the **fundamental** frequency (one cycle per period)
- `n = 2`: the **second harmonic** (two cycles per period)
- `n = 3, 4, ...`: higher harmonics

The coefficients `aₙ` and `bₙ` ARE the frequency-domain description — they tell you how much of each harmonic is hiding in the signal. They are exactly the **projections** we'll compute in Section 5.

The amazing thing: with the right `aₙ, bₙ`, you can build a square wave, a sawtooth, a heartbeat — anything periodic — out of pure sines and cosines. Add more harmonics → better approximation.

### Step 4: Enter Euler's formula — why complex exponentials?

At each frequency `n`, we track **two** real coefficients (one cosine, one sine). That's awkward: a sine and a cosine at the same frequency are really the same wave shape, just shifted by 90°. Every time you delay a wave in time, its cosine and sine amounts mix in a tangled way.

**Euler's formula** packages a cosine and a sine at the same frequency into one complex number:

```
e^(iθ) = cos(θ) + i·sin(θ)
```

Geometrically, `e^(iθ)` is a point on the **unit circle in the complex plane**, at angle θ measured from the positive real axis:

```
            Imaginary
              ▲
              │
       e^(iθ) ●
            ╱ │
          ╱   │ sin(θ)
        ╱ θ   │
      ●───────┴─→ Real
              cos(θ)
```

So:
- **Real part** of `e^(iθ)` = `cos(θ)`
- **Imaginary part** of `e^(iθ)` = `sin(θ)`

When you let θ = ωt, the point `e^(iωt)` **spins** around the unit circle as time `t` grows — at angular speed ω. Each spin contains both a cosine wave AND a sine wave at frequency ω, automatically in sync.

### Step 5: Why this packaging is so useful

Three reasons we prefer `e^(iωt)` over separate sines and cosines:

1. **Compactness.** One complex number per frequency replaces two real numbers (cos amplitude + sin amplitude). The complex magnitude is the total amplitude; the complex angle is the phase.

2. **Algebra becomes trivial.**
   - Multiplication: `e^(iA) · e^(iB) = e^(i(A+B))` — phases add.
   - Differentiation: `d/dt [e^(iωt)] = iω · e^(iωt)` — derivatives turn into multiplication by `iω`. (For sines and cosines, derivatives swap them around messily.)
   - Time-shift: shifting `t` by τ multiplies `e^(iωt)` by the constant `e^(-iωτ)`. Clean.

3. **Complex exponentials are eigenvectors of LTI systems** (remember Section 1?). For any linear time-invariant system — a filter, an amplifier, an echo chamber — `e^(iωt)` goes in and the **same** `e^(iωt)` comes out, scaled by a complex number (the frequency response at ω). That's the deep reason Fourier analysis is the natural language for signal processing.

### Step 6: From series to transform

The Fourier **series** decomposes a periodic signal into discrete harmonics (n = 0, 1, 2, ...).

But what if the signal isn't periodic — say, a single drum hit that never repeats? Then we let the period T → ∞ and the discrete harmonics squeeze together into a continuous range of frequencies. The sum becomes an integral:

```
          ∞
F(ω) = ∫ f(t)·e^(-iωt) dt
        −∞
```

This is **the Fourier transform**. It takes a (non-periodic) time-domain signal `f(t)` and gives back the frequency-domain function `F(ω)`.

Reading the formula piece by piece:
- `e^(iωt)` is the basis function for frequency ω (our "axis" in signal space).
- `e^(-iωt)` is its complex conjugate, used because the proper inner product for complex functions conjugates one factor.
- The integral is the inner product `⟨f, e^(iωt)⟩` — the **projection** of `f` onto the frequency-ω axis.

So:
> `F(ω)` = "How much of frequency ω is inside `f(t)`?"

This is the exact same idea as `c = v · e` for vectors. The Fourier transform is just projection-onto-an-axis, with the axes being complex sinusoids and the dot product being an integral.

In Section 4 we'll prove the orthogonality property that makes this all work. In Section 5 we'll compute these projections on a concrete signal.

---

## Section 4 — Why sin/cos Integrals Equal Zero

The central fact that makes Fourier analysis work:

```
∫₀²π sin(m·t) · sin(n·t) dt = 0    if m ≠ n
∫₀²π cos(m·t) · cos(n·t) dt = 0    if m ≠ n
∫₀²π sin(m·t) · cos(n·t) dt = 0    always
```

### Step 1: The intuition first (before any math)

Picture `sin(2t)` and `sin(3t)` plotted on the same graph from 0 to 2π. They oscillate at different rates. Now imagine multiplying them point by point to get a new curve `sin(2t)·sin(3t)`.

That product curve goes **positive** in some regions and **negative** in others. The integral is just the **net area** under that curve.

The key claim: when two pure sinusoids have different frequencies, the positive and negative areas of their product **cancel exactly** over a full period. The net area is zero.

When the frequencies are the **same** (m = n), the product is `sin²(mt)`, which is always ≥ 0 — no cancellation, integral is positive.

That cancellation = orthogonality.

### Step 2: Prove it with product-to-sum identities

There's a trig identity:

```
sin(A) · sin(B) = ½ [ cos(A − B) − cos(A + B) ]
```

Apply it with A = m·t, B = n·t:

```
sin(mt) · sin(nt) = ½ [ cos((m−n)t) − cos((m+n)t) ]
```

So the integral becomes:

```
∫₀²π sin(mt)·sin(nt) dt = ½ ∫₀²π cos((m−n)t) dt  −  ½ ∫₀²π cos((m+n)t) dt
```

### Step 3: Evaluate each piece

For any nonzero integer k:

```
∫₀²π cos(k·t) dt = [sin(k·t)/k]₀²π = (sin(2πk) − sin(0)) / k = 0
```

Because sin of any integer multiple of 2π is 0.

So:
- If m ≠ n, then (m−n) ≠ 0 and (m+n) ≠ 0, so **both** integrals are 0 → total is **0**. ✓
- If m = n, then (m−n) = 0, so `cos((m−n)t) = cos(0) = 1`. The first integral becomes `∫₀²π 1 dt = 2π`. The second integral is still 0. Total: ½·(2π) = **π**, not zero.

That's the full story:
```
∫₀²π sin(mt)·sin(nt) dt = { 0   if m ≠ n
                          { π   if m = n  (≠ 0)
```

The same identity machinery handles the cos·cos and sin·cos cases — they all collapse to integrals of cos(kt) or sin(kt) over a full period, which are zero.

### Step 4: Concrete numerical sanity check

Let m = 2, n = 3. The integrand is `sin(2t)·sin(3t) = ½[cos(t) − cos(5t)]`.
- `∫₀²π cos(t) dt` = sin(2π) − sin(0) = 0 ✓
- `∫₀²π cos(5t) dt` = (sin(10π) − sin(0))/5 = 0 ✓
- Sum: 0. Done.

---

## Section 5 — Projection Example: f(t) = cos(2t) + 3·sin(5t)

Let's project this signal onto several frequencies and watch orthogonality do its work.

### Setup

Our signal:
```
f(t) = cos(2t) + 3·sin(5t)
```

Our basis functions on the interval [0, 2π]:
```
cos(t), cos(2t), cos(3t), ...
sin(t), sin(2t), sin(3t), ...
```

From Section 4, we know:
```
∫₀²π cos(mt)·cos(nt) dt = π  if m = n,  else 0
∫₀²π sin(mt)·sin(nt) dt = π  if m = n,  else 0
∫₀²π sin(mt)·cos(nt) dt = 0  always
```

To extract "how much of `cos(nt)`" or "how much of `sin(nt)`" is in `f`, the projection formula is:
```
a_n = (1/π) ∫₀²π f(t) · cos(nt) dt
b_n = (1/π) ∫₀²π f(t) · sin(nt) dt
```

(The factor 1/π just normalizes — it's 1/π because the self-overlap of a sinusoid over [0, 2π] is π.)

### Projection 1 — onto cos(2t) (should hit, since signal has it with amplitude 1)

```
a_2 = (1/π) ∫₀²π [cos(2t) + 3·sin(5t)] · cos(2t) dt
```

Distribute:
```
= (1/π) ∫₀²π cos(2t)·cos(2t) dt   +   (3/π) ∫₀²π sin(5t)·cos(2t) dt
       ↑                                  ↑
   same frequency cos·cos              sin·cos (always 0)
   → π                                 → 0
```

#### What does the multiplication actually look like? (the LHS before integrating)

Before evaluating the integrals, let's see what the products `cos(2t)·cos(2t)` and `sin(5t)·cos(2t)` look like as curves. The integral is just the **net area** under each curve, so the shape tells us everything.

**Piece 1: `cos(2t) · cos(2t) = cos²(2t)`**

Use the identity `cos²(x) = (1 + cos(2x))/2`:
```
cos(2t)·cos(2t) = ½ + ½·cos(4t)
                  ↑      ↑
              DC level   pure oscillation
              (constant) (cancels to 0)
```

A constant `½` plus a wiggle. Because `cos²` is never negative, the curve stays **above the t-axis** the whole time:

```
1 ┤   ╭╮    ╭╮    ╭╮    ╭╮
  │  ╱  ╲  ╱  ╲  ╱  ╲  ╱  ╲
½ ┤───────────────────────────  ← DC level ½  (the curve never dips below 0)
  │ ╱    ╲╱    ╲╱    ╲╱    ╲
0 ┴───────────────────────────→  t
  0            π            2π
```

When we integrate, the `½·cos(4t)` part cancels (equal area above and below its own axis), but the constant `½` accumulates: `½ · 2π = π`. **Integral = π.** ✓

**Piece 2: `sin(5t) · cos(2t)`**

Use product-to-sum: `sin(A)·cos(B) = ½ [sin(A+B) + sin(A−B)]`:
```
sin(5t)·cos(2t) = ½·sin(7t) + ½·sin(3t)
                  ↑            ↑
              pure oscillation  pure oscillation
              (cancels)        (cancels)
```

**No DC term.** Both pieces are pure sines, oscillating symmetrically above and below zero:

```
  ╱╲      ╱╲╱╲      ╱╲    ╱╲      ╱╲╱╲      ╱╲
0───────────────────────────────────────────────→  t
  ╲╱      ╲╱╲╱      ╲╱    ╲╱      ╲╱╲╱      ╲╱
  0                       π                       2π
```

Equal area above and below the axis → they cancel exactly. **Integral = 0.**

**Key insight:** A non-zero integral happens **only when the product has a non-zero DC (constant) component**. That happens **only** when the two factors have matching frequency and matching shape (cos·cos or sin·sin at the same n). Every other product is pure oscillation that cancels.

So:
```
a_2 = (1/π)·π + (3/π)·0 = 1   ✓
```

We recovered the amplitude 1 exactly — coming entirely from the DC term hiding inside `cos²(2t)`.

### Projection 2 — onto sin(5t) (should hit with amplitude 3)

```
b_5 = (1/π) ∫₀²π [cos(2t) + 3·sin(5t)] · sin(5t) dt

    = (1/π) ∫₀²π cos(2t)·sin(5t) dt   +   (3/π) ∫₀²π sin(5t)·sin(5t) dt
            ↑                                       ↑
        sin·cos = 0                            same-frequency sin·sin = π

    = 0 + (3/π)·π = 3   ✓
```

Amplitude 3 recovered exactly.

### Projection 3 — onto cos(3t) (should be zero — signal has no frequency-3 cosine)

```
a_3 = (1/π) ∫₀²π [cos(2t) + 3·sin(5t)] · cos(3t) dt

    = (1/π) ∫₀²π cos(2t)·cos(3t) dt   +   (3/π) ∫₀²π sin(5t)·cos(3t) dt
            ↑                                       ↑
        cos·cos different freqs = 0            sin·cos = 0

    = 0   ✓
```

#### What do these products look like?

Both products are mismatched, so both should be pure oscillations with no DC term:

**Piece 1: `cos(2t) · cos(3t)`**

Using `cos(A)·cos(B) = ½[cos(A−B) + cos(A+B)]`:
```
cos(2t)·cos(3t) = ½·cos(t) + ½·cos(5t)
                  ↑           ↑
              pure oscillation  pure oscillation
              (cancels)        (cancels)
```

No DC. Picture two clean cosines (at frequencies 1 and 5) added together — wiggles symmetric around zero. Integral cancels to 0.

**Piece 2: `sin(5t) · cos(3t)`**

Using `sin(A)·cos(B) = ½[sin(A+B) + sin(A−B)]`:
```
sin(5t)·cos(3t) = ½·sin(8t) + ½·sin(2t)
                  ↑           ↑
              pure oscillation  pure oscillation
              (cancels)        (cancels)
```

Again — two pure sines added together, no DC term, integral cancels to 0.

```
  Mismatched product (e.g. cos(2t)·cos(3t)):

    ╱╲╱╲    ╱╲      ╱╲╱╲    ╱╲╱╲      ╱╲    ╱╲╱╲
  0─────────────────────────────────────────────→  t
    ╲╱╲╱    ╲╱      ╲╱╲╱    ╲╱╲╱      ╲╱    ╲╱╲╯
    0                       π                       2π

  Equal positive and negative area → net integral = 0.
```

Compare to Projection 1's `cos²(2t)` curve, which stayed entirely **above** the t-axis (because of the hidden DC term `½`). That's the visual difference between a matching projection and a mismatched one.

### Projection 4 — onto sin(2t) (should be zero — signal has cos(2t), not sin(2t))

This one is interesting. The signal has frequency 2 in it, but only as a **cosine**, not a sine. Cos and sin of the **same** frequency are still orthogonal:

```
b_2 = (1/π) ∫₀²π [cos(2t) + 3·sin(5t)] · sin(2t) dt

    = (1/π) ∫₀²π cos(2t)·sin(2t) dt   +   (3/π) ∫₀²π sin(5t)·sin(2t) dt
            ↑                                       ↑
        sin·cos = 0                           different freq sin·sin = 0

    = 0   ✓
```

Zero again — orthogonality keeps cos and sin components separate even at the same frequency.

### One integral, fully cranked by hand

We claimed `∫₀²π cos(2t)·cos(3t) dt = 0`. Let's actually evaluate it.

Use the product-to-sum identity:
```
cos(A)·cos(B) = ½ [cos(A − B) + cos(A + B)]
```

With A = 2t, B = 3t:
```
cos(2t)·cos(3t) = ½ [cos(−t) + cos(5t)] = ½ [cos(t) + cos(5t)]
```
(since cos is even)

Integrate:
```
∫₀²π cos(2t)·cos(3t) dt = ½ ∫₀²π cos(t) dt  +  ½ ∫₀²π cos(5t) dt

                       = ½ [sin(t)]₀²π    +  ½ [sin(5t)/5]₀²π

                       = ½ · (sin(2π) − sin(0))   +   ½ · (sin(10π) − sin(0))/5

                       = ½ · (0 − 0)              +   (1/10) · (0 − 0)

                       = 0   ✓
```

Cancellation happens because we integrate full periods of cos(t) and cos(5t) — the positive and negative halves wipe out.

And the matching case `∫₀²π cos²(2t) dt`:
```
cos²(2t) = (1 + cos(4t))/2

∫₀²π cos²(2t) dt = ½ ∫₀²π 1 dt   +   ½ ∫₀²π cos(4t) dt

                = ½ · 2π          +   ½ · 0

                = π   ✓
```

The constant part `½` has nothing to cancel — it accumulates to `π`. That's why matching-frequency projections survive while everything else dies.

### Full Fourier picture of f(t)

If we project `f(t) = cos(2t) + 3·sin(5t)` onto **every** basis frequency `n = 0, 1, 2, 3, ...`:

| n | a_n (cos coeff) | b_n (sin coeff) |
|---|-----------------|-----------------|
| 0 | 0               | —               |
| 1 | 0               | 0               |
| 2 | **1**           | 0               |
| 3 | 0               | 0               |
| 4 | 0               | 0               |
| 5 | 0               | **3**           |
| 6 | 0               | 0               |
| ... | 0             | 0               |

**Every** coefficient is zero except the two that match the signal's actual content. That's the Fourier transform in its purest form: a sparse list of "how much of each pure tone is in f."

And from that table you could reconstruct `f(t)` exactly:
```
f(t) = 1·cos(2t) + 3·sin(5t)
```
which is what we started with.

### What this demonstrates

1. **Matching frequency + matching shape (cos/cos or sin/sin)** → integral = π → you recover the amplitude.
2. **Any mismatch** (different frequency, OR cos vs sin at same frequency) → integral = 0 → no leakage.

This is why the Fourier transform is the "right" coordinate system for signals: each frequency component is read off independently, with zero interference from the others. That's orthogonality earning its keep.

You'll see exactly the table above — zeros everywhere except `a_2 = 1` and `b_5 = 3`.
