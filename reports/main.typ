/////////////
// Imports //
/////////////
#import "@preview/adaptable-pset:0.2.0": *
#import "@preview/physica:0.9.8": *
#import "@preview/unify:0.8.1": *
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#show: codly-init.with()
#codly(languages: codly-languages)

/////////////////
// Maths Setup //
/////////////////

// upright vectors
#let vectorboldupright(a) = vb($upright(#a)$)
#let vbu = vectorboldupright
#let vectorunitupright(a) = vu($upright(#a)$)
#let vuu = vectorunitupright
#let vectorarrowupright(a) = va($upright(#a)$)
#let vau = vectorarrowupright

// automatically use square brackets for vectors and matricies
#set math.vec(delim: "[")
#set math.mat(delim: "[")
#let vecrowOld = vecrow
#let vecrow = vecrowOld.with(delim: "[")

////////////////////
// Document Setup //
////////////////////

// assignment info
#show: homework.with(
    title: "HW02",
    author: "Vai Srivastava",
    collaborators: [],
    course-id: "ENAE 641: Linear System Dynamics",
    instructor: "Dr. Robert Sanner",
    semester: "Fall 2026",
    due-time: datetime(
        // due date
        year: 2026,
        month: 10,
        day: 13,

        // due time
        hour: 23,
        minute: 59,
        second: 00,
    ).display("[month repr:long] [day padding:zero], [year] at [hour repr:24]:[minute padding:zero]:[second padding:zero]"),

    // (defaults to A4)
    paper-size: "us-letter",
)

// document settings
#set text(font: "New Computer Modern", size: 10pt)
#set enum(numbering: "a)")

// problem headings
#let probOld = prob
#let prob = prob.with(color: black)

////////////////////////////
// The Assignment Itself: //
// Problems and Solutions //
////////////////////////////

#prob[
    For a time-invariant linear system with:
    $
        vb(A) = mat(
            0, 0, -13;
            1, 0, -17;
            0, 1, -5;
        ) quad
        vb(B) = mat(
            13, 1;
            4, 1;
            1, 0;
        ) quad
        vb(C) = mat(
            1, -1, 1;
            0, 1, -1;
        ) quad
    $
    and $vb(D) = vb(0)$:

    + Determine the representation of the state dynamics in the real, block-modal form, and identify the required coordinate transformation matrix $vb(P)$.
    <hwk:p01a>

    + Give a physically meaningful (i.e. real valued) initial state $vbu(x)_0$ which will excite #strong[only] the oscillatory mode of the homogeneous state dynamics. Repeat for the purely exponential mode.
    <hwk:p01b>

    + Explicitly compute the transfer function matrix $vb(G)(s)$ using the transformed state-space model in #link(<hwk:p01a>)[Part A], and verify that the result is identical to that obtained using the original state-space model.
    <hwk:p01c>
    // Python can be used for this problem
] <hwk:p01>

// TODO: answer
+ Answer
<hwk:s01a>

// TODO: answer
+ Answer
<hwk:s01b>

// TODO: answer
+ Answer
<hwk:s01c>

#pagebreak(weak: true)

#prob[
    The homogeneous state equation for a linear system has matrix $vb(A)(t)$ that can be partitioned as
    $
        mat(
            vb(A)_1(t), vb(A)_(1, 2)(t);
            0, vb(A)_2(t);
        )
    $
    $vb(A)_1$ and $vb(A)_2$ are square, but of possibly different sizes, while $vb(A)_(1,2)$ is not necessarily square.

    + Show that the state transition matrix has the compatible form:
        $
            vb(Phi)(t, tau) = mat(
                vb(Phi)_1(t), vb(Phi)_(1, 2)(t);
                0, vb(Phi)_2(t);
            )
        $
        where $vb(Phi)_1$ is the STM for $vb(A)_1$ and $vb(Phi)_2$ s the STM for $vb(A)_2$
    <hwk:p02a>

    + Give an expression for $vb(Phi)_(1, 2)(t, tau)$ in terms of $vb(Phi)_1$, $vb(Phi)_2$, and $vb(A)_(1, 2)$.

        #smallcaps[Note:] $vb(Phi)_(1, 2)$ will not be a STM for $vb(A)_(1, 2)$, since that matrix is not necessarily square.
    <hwk:p02b>
] <hwk:p02>

// TODO: answer
+ Answer
<hwk:s02a>

// TODO: answer
+ Answer
<hwk:s02b>

#pagebreak(weak: true)

#prob[
    Given the state-space model:
    $
        dot(x)_1(t) &= t x_2(t) \
        dot(x)_2(t) &= -x_2(t) + e^(-t) u(t) \
        y(t) &= t x_1(t)
    $

    + Find the state transition matrix $Phi(t, tau)$.
    <hwk:p03a>

    + Find the impulse response $G(t, tau)$.
    <hwk:p03b>

    + Find the response $y(t)$ for $t gt.eq 1$ when $vb(x)(1)^TT = vecrow(1, -2)$ and $u(t) = 1$ for all $t gt.eq 1$.
    <hwk:p03a>
] <hwk:p03>

// TODO: answer
+ Answer
<hwk:s03a>

// TODO: answer
+ Answer
<hwk:s03b>

// TODO: answer
+ Answer
<hwk:s03c>

#pagebreak(weak: true)

#prob[
    The characteristic polynomial of the $vb(A)$ matrix for an LTI dynamic system can be written as $r(lambda) = (lambda^2 + 2 lambda + 10)(lambda + 3)^3$. You also know that the rank of $(vb(A) + 3 vb(I))$ is 4.

    Let $vb(P)$ be the matrix which will transform $vb(A)$ to an appropriate real block modal form. Assume the $vb(P)$ matrix is such that the components of the modal form arising from the complex eigenvalues will be located in the upper left corner.

    + Compute the matrix $vb(P)^(-1) e^(vb(A) t) vb(P)$.
    <hwk:p04a>

    + Suppose that we identify the outputs of this system to be $vb(y)(t) = vb(P)^(-1) vbu(x)(t)$. If the initial condition is $vbu(x)_0 = vbu(p)_4$ where $vbu(p)_4$ is the fourth column of the $vb(P)$ matrix, find $vbu(y)(t)$ explicitly as a function of time.
    <hwk:p04b>

    + Suppose in #link(<hwk:p04b>)[Part B] that the inital condition is instead $vbu(x)_0 = vbu(p)_2$. Find $vbu(y)(t)$.
    <hwk:p04c>
] <hwk:p04>

// TODO: answer
+ Answer
<hwk:s04a>

// TODO: answer
+ Answer
<hwk:s04b>

// TODO: answer
+ Answer
<hwk:s04c>

#pagebreak(weak: true)

#prob[
    Suppose the matrix $vb(A)(t)$ for a linear system can be written as
    $
        vb(A)(t) = sum_(j = 1)^r f_j (t) vb(A)_j
    $
    where the $vb(A)_j$ are constant, and additionally $vb(A)_i vb(A)_j = vb(A)_j vb(A)_i$ (commutative) for every $i$, $j$.

    + Show that the state transition matrix for this system can be written as a matrix exponential function.
    <hwk:p05a>

    + Show that the matrix exponential form for the STM can in fact be factored as
        $
            vb(Phi)(t, tau) = e^(vb(A)_1) integral_tau^t f_1(sigma) dd(sigma) e^(vb(A)_2) integral_tau^t f_2(sigma) dd(sigma) space ... space e^(vb(A)_r) integral_tau^t f_r(sigma) dd(sigma)
        $
    <hwk:p05b>

    + Use the observations above to compute the state transition matrix for
        $
            vb(A)(t) = mat(
                cos(3t), sin(3t);
                -sin(3t), cos(3t);
            )
        $
    <hwk:p05c>
] <hwk:p05>

// TODO: answer
+ Answer
<hwk:s05a>

// TODO: answer
+ Answer
<hwk:s05b>

// TODO: answer
+ Answer
<hwk:s05c>

#pagebreak(weak: true)

#prob[
    For the homogeneous linear system $dot(vbu(x))(t) = vb(A)(t) vbu(x)(t)$, find the state transition matrix $vb(Phi)(t, tau)$ when:

    +
        $
            vb(A)(t) = mat(
                t, 0;
                0, cos(t);
            )
        $
    <hwk:p06a>

    +
        $
            vb(A)(t) = mat(
                cos(2t), cos(2t);
                0, 0;
            )
        $
    <hwk:p06b>

    +
        $
            vb(A)(t) = mat(
                -2, sin(t);
                0, -1;
            )
        $
    <hwk:p06c>

    +
        $
            vb(A)(t) = mat(
                -3, -1;
                -1, -3;
            )
        $
    <hwk:p06d>

    Show all work for your solutions. In #link(<hwk:p06d>)[Part D], use a dyadic expansion approach.
] <hwk:p06>

// TODO: answer
+ Answer
<hwk:s06a>

// TODO: answer
+ Answer
<hwk:s06b>

// TODO: answer
+ Answer
<hwk:s06c>

// TODO: answer
+ Answer
<hwk:s06d>

#pagebreak(weak: true)

#prob[
    Consider a time-periodic system with
    $
        vb(A)(t) = mat(
            -1 + 1.5 cos^2(t), 1 - 1.5 sin(t) cos(t);
            -1 - 1.5 sin(t) cos(t), -1 + 1.5 sin^2(t);
        )
    $

    + Find the eigenvalues of $vb(A)(t)$, show that they do not change with time, and are located in the left-half of the complex plane.
    <hwk:p07a>

    + Verify that the two fundamental solutions for this system are:
        $
            vb(psi)_1(t) = mat(
                e^(-t) sin(t);
                e^(-t) cos(t);
            ) quad
            vb(psi)_2(t) = mat(
                e^(t/2) cos(t);
                -e^(t/2) sin(t);
            )
        $
    <hwk:p07b>

    + Determine the state transition matrix for this system. Is this system stable; i.e. does $vbu(x)(t) -> vb(0)$ from every initial condition $vbu(x)_0$?
    <hwk:p07c>

    + Comparing your results in #link(<hwk:p07a>)[Part A] and #link(<hwk:p07c>)[Part C], comment on the utility of the eigenvalues of a time-varying $vb(A)(t)$ for analyzing stability.
    <hwk:p07d>
] <hwk:p07>

// TODO: answer
+ Answer
<hwk:s07a>

// TODO: answer
+ Answer
<hwk:s07b>

// TODO: answer
+ Answer
<hwk:s07c>

// TODO: answer
+ Answer
<hwk:s07d>

#pagebreak(weak: true)

== Code

#codly(header: [./src/index.py])
#raw(read("../src/index.py"), block: true, lang: "python") <code:index.py>
