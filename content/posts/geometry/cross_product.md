---
title: "Useful Properties of Cross Product"
date: 2020-12-27T12:55:13-08:00
draft: false
tags: ["geometry", "navigation"]
---
This note is going to go important properties of what you might call the cross product matrix.
For a vector $v \in \mathbb{R}^3$ we'll write $$[v]_\times = \begin{bmatrix} 0 & -v_z & v_y \\\\ v_z & 0 & -v_x \\\\ -v_y & v_x & 0\end{bmatrix}$$ to mean the $3\times 3$ matrix such that $[v]\_\times w = v \times w$ for a all $w \in \mathbb{R}^3$.

These matrices are closely related to rotations due through the Rodrigues formula:
$$
R(\theta, u) = \exp\left([\theta u ]\_\times\right) = I + \sin(\theta) [u]\_\times + (1 - \cos(\theta))[u]\_\times^2.
$$
Here the matrix $R(\theta, a)$ is a rotation about the unit norm vector $a$ by an angle $\theta$.
# Relationship to Linear Transforms
Important properties of the cross-product matrix can be derived from its relationship two the determinant:
$$
\det\left(\begin{bmatrix}
u & v & w
\end{bmatrix}\right) = u^T [v]_\times w
\qquad \forall \; u, v, w \in \mathbb{R}^3.$$
This can be easily proven using expansion by minors along the first column.
Given an arbitrary invertible matrix $M \in \mathbb{R}^{3 \times 3}$ we therefore have:

$$
\det\left(M\begin{bmatrix}
u & v & w
\end{bmatrix}\right) = u^T\left(M^T [Mv]_\times M\right)w
\qquad \forall \; u, v, w \in \mathbb{R}^3.$$

But basic properties of the determinant also give use
$$
\det\left(M\begin{bmatrix}
u & v & w
\end{bmatrix}\right) = \det(M) u^T [v]_\times w
\qquad \forall \; u, v, w \in \mathbb{R}^3.$$
Since this holds for all $u$ and $v$ we can conclude that:

$$ \det(M)  M^{-T} [v]_\times = [ Mv ]\_\times M \qquad \forall\; v \in \mathbb{R}^3, M \in GL(3).$$

This identity is particularly useful when $M = R$ where $R$ is an element of $SO(3)$ so that $R^T = R^{-1}$ and $\det(R) = 1$:

$$R [v]_\times = [ Rv ]\_\times R \qquad \forall \; v \in \mathbb{R}^3, R \in SO(3).$$

# Algebraic Properties
## Anti-symmetry
It is easy to verify that:
$$
[v]\_\times u = -[u]\_\times v.
$$

## Lie Bracket
The following relationship is frequently applied in navigation equations:
$$ [[u]\_\times v]\_\times = vu^T - uv^T $$
It can be verified explicitly (as a hint, the anti-symmetry and bilinearity of the expression means it is sufficient to test the result for the pairs $(u, v) \in \\{(e_x, e_y), (e_y, e_z), (e_x, e_z)\\}$).

## Power Relationships
One of the most useful identities working with cross-product matrices is:
$$
\begin{equation}
[u]_\times [v]\_\times = vu^T - u^T v I
\end{equation}
$$
This is not too painful to verify explicitly.  We see from this and the Lie Bracket property that:

$$
\begin{equation}
[u]_\times [v]\_\times - [v]\_\times [u]\_\times = vu^T - uv^T =\left[[u]\_\times v\right]\_\times
\end{equation}
$$

For a unit vector $u \in \mathbb{R}^3$ Equation (1) establishes an important algebraic relationship:

$$
\begin{align}
[u ]_\times^3 &= (uu^T - I) [u]\_\times = -[u]\_\times \\\\
[u ]\_\times^4 &= -[u]\_\times^3 [u]\_\times = -[u]\_\times^2
\end{align}
$$
This makes it clear that any power series in $[u]_x$ can always be rewritten in at most three terms:

$$
\alpha I + \beta [u]\_\times + \gamma [u]\_\times^2
$$

I find it more useful to write such expressions as:

$$
a I + b [u]\_\times + c uu^T
$$
with $a = \alpha - \gamma$, $b = \beta$ and $c = \gamma$.

It is also useful to note that the inverse of such a matrix (when it exists) can be written:
$$
(a I + b [u]\_\times + c uu^T)^{-1} = d I + e [u]\_\times + f uu^T,
$$
with:
$$
\begin{bmatrix}d \\\\ e \\\\ f\end{bmatrix}
= \frac{1}{a^2+b^2}\begin{bmatrix}a \\\\ -b \\\\ \frac{b^2-ca}{a+c}\end{bmatrix}.
$$

> *Proof:* We explicitly expand the product of these two matrices:
$$
\begin{gather}
(a I + b [u]\_\times + cuu^T) (d I + e [u]\_\times + f uu^T) \\\\
=\\\\
ad I + ae [u]\_\times + af uu^T + bd [u]\_\times + be uu^T - be I + c(d+f) uu^T
\end{gather}
$$
For the second matrix to be the inverse of the first, we require that:
$$
\begin{bmatrix}
a & -b & 0 \\\\
b &  a& 0 \\\\
c & b & a+c
\end{bmatrix}
\begin{bmatrix}
d \\\\ e \\\\ f
\end{bmatrix} = \begin{bmatrix}1 \\\ 0 \\\ 0\end{bmatrix}
$$
The result can then be found by using the general matrix identity:
$$
\begin{bmatrix}
A & 0 \\\\
B & C
\end{bmatrix}^{-1} = \begin{bmatrix}
A^{-1} & 0 \\\\
-C^{-1}BA^{-1} & C^{-1}
\end{bmatrix}.
$$



# Eigenvector Analysis
First, it is easy to verify that $[u]_\times u = 0$ so that $u$
itself is an eigenvector with eigenvalue 0.

If $u = 0$ then all the eigenvalues are trivially zero.
Otherwise,  let $v$ be any vector orthogonal to $u$ and let $w =\frac{1}{\\|u\\|} [u]_\times v$.
Then $u + j w$ is an eigenvector with eigenvalue $-j\\|v\\|$:
$$
\begin{align}
[u]_x (v + jw) &= \\|u\\| w + j [u]_x w\\\\
&=  \\|u\\|w + j \frac{1}{\\|u\\|}[u]_x^2 v \\\\
&=  \\|u\\|w + j \frac{1}{\\|u\\|}(uu^T - \\|u\\|^2 I) v \\\\
&=  \\|u\\|\left(w - j  v \right) = -j\\|u\\|\left(v+jw\right).
\end{align}
$$
That $v - jw$ is an eigenvector with eigenvalue $-j\\|u\\|$ can be shown similarly.

This line or reasoning can be made more concrete by analyzing the
matrix $[e_x]\_\times$ with $u = e_y$ and $w = e_z$, and then
concluding that the resulting analysis can be transfered to an
arbitrary $[v]\_\times$ by a rotation matrix and scaling.

Armed with this eigenvector analysis we can revisit the matrices of the form:
$$
aI + b [u]\_\times + cuu^T.
$$
We see that these matrices have the same eigenvectors as $[u]\_\times$ owing to $u^Tv = u^Tw = 0$.
When $\|u\| = 1$ we can see that the eigen-vector / eigen-value pairs are:
- $(u, a + c)$
- $(v + jw, a-bj)$
- $(v - jw, a+bj)$

This leads to a second approach to finding our general equation for
the inverse of these matrices, based on inverting these eigenvalues
and term-matching to find $(d, e, f)$ from the prior section.
