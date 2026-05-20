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

---

## Section 3 — Why sin/cos Integrals Equal Zero

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

You can verify in Python:
```python
from scipy.integrate import quad
from math import sin, pi
quad(lambda t: sin(2*t)*sin(3*t), 0, 2*pi)  # ≈ 0
quad(lambda t: sin(2*t)*sin(2*t), 0, 2*pi)  # ≈ π
```

---

## Section 4 — Projection Example: f(t) = cos(2t) + 3·sin(5t)

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

From Section 3, we know:
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

So:
```
a_2 = (1/π)·π + (3/π)·0 = 1   ✓
```

We recovered the amplitude 1 exactly.

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

Cleanly zero. The signal has no frequency-3 cosine content.

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

### Numerical verification

```python
from scipy.integrate import quad
from math import cos, sin, pi

f = lambda t: cos(2*t) + 3*sin(5*t)

for n in range(7):
    a_n, _ = quad(lambda t: f(t)*cos(n*t), 0, 2*pi)
    b_n, _ = quad(lambda t: f(t)*sin(n*t), 0, 2*pi)
    print(f"n={n}: a_n={a_n/pi:.4f}, b_n={b_n/pi:.4f}")
```

You'll see exactly the table above — zeros everywhere except `a_2 = 1` and `b_5 = 3`.
