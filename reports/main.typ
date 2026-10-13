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
    ).display(
        "[month repr:long] [day padding:zero], [year] at [hour repr:24]:[minute padding:zero]:[second padding:zero]",
    ),

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

+ Real block-modal form

    $vb(A)$ is a companion matrix, so its characteristic polynomial can be read off its last column:
    $
        det(s vb(I) - vb(A)) = s^3 + 5 s^2 + 17 s + 13 = (s + 1)(s^2 + 4 s + 13)
    $
    The eigenvalues are $lambda_(1,2) = -2 plus.minus 3j$ (the oscillatory mode) and $lambda_3 = -1$ (the purely exponential mode).

    For an eigenvalue $lambda$, the rows of $(vb(A) - lambda vb(I)) vbu(v) = vb(0)$ are
    $
        -lambda v_1 - 13 v_3 = 0, quad
        v_1 - lambda v_2 - 17 v_3 = 0, quad
        v_2 - (lambda + 5) v_3 = 0
    $
    Setting $v_3 = 1$, the last two rows give $v_2 = lambda + 5$ and $v_1 = lambda^2 + 5 lambda + 17$. The first row then reads $-(lambda^3 + 5 lambda^2 + 17 lambda + 13) = 0$, which holds because $lambda$ is a root of the characteristic polynomial. Thus
    $
        vbu(v)(lambda) = vec(lambda^2 + 5 lambda + 17, lambda + 5, 1)
        quad => quad
        vbu(v)_1 = vbu(v)(-2 + 3j) = vec(2, 3, 1) + j vec(3, 3, 0), quad
        vbu(v)_3 = vbu(v)(-1) = vec(13, 4, 1)
    $

    Write $vbu(v)_1 = vbu(a) + j vbu(b)$ and $lambda_1 = sigma + j omega$ with $sigma = -2$, $omega = 3$. The real and imaginary parts of $vb(A) vbu(v)_1 = lambda_1 vbu(v)_1$ are
    $
        vb(A) vbu(a) = sigma vbu(a) - omega vbu(b), quad
        vb(A) vbu(b) = omega vbu(a) + sigma vbu(b)
        quad => quad
        vb(A) mat(vbu(a), vbu(b)) = mat(vbu(a), vbu(b)) mat(sigma, omega; -omega, sigma)
    $
    so $vbu(a)$ and $vbu(b)$ span a real invariant subspace carrying the oscillatory mode. Placing the complex block in the upper left, the coordinate transformation $vbu(x) = vb(P) vbu(z)$ uses
    $
        vb(P) = mat(Re(vbu(v)_1), Im(vbu(v)_1), vbu(v)_3) = mat(2, 3, 13; 3, 3, 4; 1, 0, 1) quad qed
    $
    with $det(vb(P)) = -30$ and
    $
        vb(P)^(-1) = 1/30 mat(-3, 3, 27; -1, 11, -31; 3, -3, 3)
    $
    In the new coordinates, $dot(vbu(z)) = vb(Lambda) vbu(z) + macron(vb(B)) vbu(u)$ and $vbu(y) = macron(vb(C)) vbu(z)$, with $macron(vb(D)) = vb(D) = vb(0)$ and
    $
        vb(Lambda) = vb(P)^(-1) vb(A) vb(P) = mat(-2, 3, 0; -3, -2, 0; 0, 0, -1) quad qed
    $
    $
        macron(vb(B)) = vb(P)^(-1) vb(B) = mat(0, 0; 0, 1/3; 1, 0), quad
        macron(vb(C)) = vb(C) vb(P) = mat(0, 0, 10; 2, 3, 3) quad qed
    $
    The Python calculation confirms these matrices:
    #raw(read("../outputs/text/s01a.txt"), block: true)
<hwk:s01a>

+ Single-mode initial states

    The homogeneous response is $vbu(x)(t) = e^(vb(A) t) vbu(x)_0 = vb(P) e^(vb(Lambda) t) vbu(z)_0$ with $vbu(z)_0 = vb(P)^(-1) vbu(x)_0$, where
    $
        e^(vb(Lambda) t) = mat(
            e^(-2t) cos(3t), e^(-2t) sin(3t), 0;
            -e^(-2t) sin(3t), e^(-2t) cos(3t), 0;
            0, 0, e^(-t);
        )
    $
    Since $e^(vb(Lambda) t)$ is block diagonal, a mode is excited only if $vbu(z)_0$ has a nonzero component in that mode's block. Any real, nonzero $vbu(x)_0 = alpha vbu(p)_1 + beta vbu(p)_2$ has $z_3(0) = 0$ and excites only the oscillatory mode. Choosing $vbu(x)_0 = vbu(p)_1$ (so $vbu(z)_0 = vecrow(1, 0, 0)^TT$):
    $
        vbu(x)(t) = e^(-2t) (cos(3t) vbu(p)_1 - sin(3t) vbu(p)_2)
        = e^(-2t) vec(2 cos(3t) - 3 sin(3t), 3 cos(3t) - 3 sin(3t), cos(3t))
    $
    $
        vbu(x)_0 = vec(2, 3, 1) quad qed
    $
    Similarly, any nonzero real multiple of $vbu(p)_3$ has $z_1(0) = z_2(0) = 0$ and excites only the exponential mode. Choosing $vbu(x)_0 = vbu(p)_3$ (so $vbu(z)_0 = vecrow(0, 0, 1)^TT$):
    $
        vbu(x)(t) = e^(-t) vec(13, 4, 1), quad
        vbu(x)_0 = vec(13, 4, 1) quad qed
    $
    Simulating both homogeneous responses confirms that the unexcited modal coordinates stay at zero:
    #figure(
        image("../outputs/figures/s01b.png", width: 100%),
        caption: [Homogeneous response (top) and modal coordinates $vbu(z) = vb(P)^(-1) vbu(x)$ (bottom) for each single-mode initial state.],
    )
    #raw(read("../outputs/text/s01b.txt"), block: true)
<hwk:s01b>

+ Transfer function matrix

    Since $macron(vb(D)) = vb(0)$, the transformed model gives $vb(G)(s) = macron(vb(C)) (s vb(I) - vb(Lambda))^(-1) macron(vb(B))$. The matrix $s vb(I) - vb(Lambda)$ is block diagonal, so it can be inverted block by block:
    $
        mat(s + 2, -3; 3, s + 2)^(-1) = 1/(s^2 + 4s + 13) mat(s + 2, 3; -3, s + 2), quad
        (s + 1)^(-1) = 1/(s + 1)
    $
    Then
    $
        (s vb(I) - vb(Lambda))^(-1) macron(vb(B))
        = mat(
            (s + 2)/(s^2 + 4s + 13), 3/(s^2 + 4s + 13), 0;
            -3/(s^2 + 4s + 13), (s + 2)/(s^2 + 4s + 13), 0;
            0, 0, 1/(s + 1);
        ) mat(0, 0; 0, 1/3; 1, 0)
        = mat(
            0, 1/(s^2 + 4s + 13);
            0, (s + 2)/(3(s^2 + 4s + 13));
            1/(s + 1), 0;
        )
    $
    and premultiplying by $macron(vb(C))$:
    $
        macron(vb(C)) (s vb(I) - vb(Lambda))^(-1) macron(vb(B))
        = mat(0, 0, 10; 2, 3, 3) mat(
            0, 1/(s^2 + 4s + 13);
            0, (s + 2)/(3(s^2 + 4s + 13));
            1/(s + 1), 0;
        )
        = mat(
            10/(s + 1), 0;
            3/(s + 1), (2 + (s + 2))/(s^2 + 4s + 13);
        )
    $
    $
        vb(G)(s) = mat(
            10/(s + 1), 0;
            3/(s + 1), (s + 4)/(s^2 + 4s + 13);
        ) quad qed
    $
    The modal form makes the structure visible. $macron(vb(B))$ routes $u_1$ only into the exponential mode $z_3$ and $u_2$ only into the oscillatory mode, and $y_1$ observes only $z_3$ (the first row of $macron(vb(C))$ is $vecrow(0, 0, 10)$). Hence $G_(12)(s) = 0$, and the oscillatory poles do not appear in the first column.

    With the original model, every entry of $vb(C)(s vb(I) - vb(A))^(-1) vb(B)$ has the common denominator $det(s vb(I) - vb(A)) = (s + 1)(s^2 + 4s + 13)$. Using `ss` and `tf` from python-control (the Python counterparts of the Matlab commands) gives the numerators below; factoring them:
    $
        G_(11)(s) & = (10 s^2 + 40 s + 130)/((s + 1)(s^2 + 4s + 13)) = (10(s^2 + 4s + 13))/((s + 1)(s^2 + 4s + 13)) = 10/(s + 1) \
        G_(21)(s) & = (3 s^2 + 12 s + 39)/((s + 1)(s^2 + 4s + 13)) = (3(s^2 + 4s + 13))/((s + 1)(s^2 + 4s + 13)) = 3/(s + 1) \
        G_(22)(s) & = (s^2 + 5 s + 4)/((s + 1)(s^2 + 4s + 13)) = ((s + 1)(s + 4))/((s + 1)(s^2 + 4s + 13)) = (s + 4)/(s^2 + 4s + 13)
    $
    and $G_(12)(s) = 0$. After the pole-zero cancellations, this is identical to the modal-form result. This is expected, since a transfer function is invariant under a similarity transformation:
    $
        macron(vb(C)) (s vb(I) - vb(Lambda))^(-1) macron(vb(B))
        = vb(C) vb(P) (vb(P)^(-1) (s vb(I) - vb(A)) vb(P))^(-1) vb(P)^(-1) vb(B)
        = vb(C) (s vb(I) - vb(A))^(-1) vb(B)
    $
    The Python calculation, which also compares both models and the hand result along $s = j omega$:
    #raw(read("../outputs/text/s01c.txt"), block: true)
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
                vb(Phi)_1(t, tau), vb(Phi)_(1, 2)(t, tau);
                0, vb(Phi)_2(t, tau);
            )
        $
        where $vb(Phi)_1$ is the STM for $vb(A)_1(t)$ and $vb(Phi)_2$ is the STM for $vb(A)_2(t)$
    <hwk:p02a>

    + Give an expression for $vb(Phi)_(1, 2)(t, tau)$ in terms of $vb(Phi)_1$, $vb(Phi)_2$, and $vb(A)_(1, 2)$.

        #smallcaps[Note:] $vb(Phi)_(1, 2)$ will not be a STM for $vb(A)_(1, 2)$, since that matrix is not necessarily square.
    <hwk:p02b>
] <hwk:p02>

+ Block-triangular form of the STM

    The STM is the unique solution of
    $
        pdv(, t) vb(Phi)(t, tau) = vb(A)(t) vb(Phi)(t, tau), quad vb(Phi)(tau, tau) = vb(I)
    $
    Partition $vb(Phi)$ conformally with $vb(A)$, writing its blocks as $vb(X)_(i j)(t, tau)$, so that $vb(X)_(11)$ has the size of $vb(A)_1$ and $vb(X)_(22)$ has the size of $vb(A)_2$:
    $
        pdv(, t) mat(vb(X)_(11), vb(X)_(12); vb(X)_(21), vb(X)_(22))
        = mat(vb(A)_1, vb(A)_(1,2); vb(0), vb(A)_2) mat(vb(X)_(11), vb(X)_(12); vb(X)_(21), vb(X)_(22)),
        quad
        mat(vb(X)_(11), vb(X)_(12); vb(X)_(21), vb(X)_(22)) (tau, tau) = mat(vb(I), vb(0); vb(0), vb(I))
    $
    Multiplying out the blocks gives four equations:
    $
        pdv(vb(X)_(21), t) & = vb(A)_2(t) vb(X)_(21), & quad vb(X)_(21)(tau, tau) & = vb(0) \
        pdv(vb(X)_(22), t) & = vb(A)_2(t) vb(X)_(22), & quad vb(X)_(22)(tau, tau) & = vb(I) \
        pdv(vb(X)_(11), t) & = vb(A)_1(t) vb(X)_(11) + vb(A)_(1,2)(t) vb(X)_(21), & quad vb(X)_(11)(tau, tau) & = vb(I) \
        pdv(vb(X)_(12), t) & = vb(A)_1(t) vb(X)_(12) + vb(A)_(1,2)(t) vb(X)_(22), & quad vb(X)_(12)(tau, tau) & = vb(0)
    $
    The first two equations do not involve the top blocks:
    - The equation for $vb(X)_(21)$ is linear and homogeneous with a zero initial value. By uniqueness, $vb(X)_(21)(t, tau) = vb(Phi)_2(t, tau) vb(X)_(21)(tau, tau) = vb(0)$ for all $t$.
    - The equation for $vb(X)_(22)$ is exactly the defining problem for the STM of $vb(A)_2(t)$, so $vb(X)_(22) = vb(Phi)_2(t, tau)$.
    - With $vb(X)_(21) = vb(0)$, the equation for $vb(X)_(11)$ reduces to $pdv(vb(X)_(11), t) = vb(A)_1(t) vb(X)_(11)$ with $vb(X)_(11)(tau, tau) = vb(I)$. This is the defining problem for the STM of $vb(A)_1(t)$, so $vb(X)_(11) = vb(Phi)_1(t, tau)$.

    Naming the remaining block $vb(Phi)_(1,2) := vb(X)_(12)$:
    $
        vb(Phi)(t, tau) = mat(vb(Phi)_1(t, tau), vb(Phi)_(1,2)(t, tau); vb(0), vb(Phi)_2(t, tau)) quad qed
    $
<hwk:s02a>

+ Off-diagonal block

    With $vb(X)_(22) = vb(Phi)_2$, the last block equation is a forced linear equation for $vb(Phi)_(1,2)$. Its dynamics are governed by $vb(A)_1$, and its "input" is $vb(A)_(1,2)(t) vb(Phi)_2(t, tau)$:
    $
        pdv(vb(Phi)_(1,2), t) = vb(A)_1(t) vb(Phi)_(1,2) + vb(A)_(1,2)(t) vb(Phi)_2(t, tau), quad
        vb(Phi)_(1,2)(tau, tau) = vb(0)
    $
    Apply the variation-of-constants formula column by column. The homogeneous term $vb(Phi)_1(t, tau) vb(Phi)_(1,2)(tau, tau)$ vanishes, leaving
    $
        vb(Phi)_(1,2)(t, tau) = integral_tau^t vb(Phi)_1(t, sigma) vb(A)_(1,2)(sigma) vb(Phi)_2(sigma, tau) dd(sigma) quad qed
    $
    Check: at $t = tau$ the integral vanishes. By the Leibniz rule, with $pdv(, t) vb(Phi)_1(t, sigma) = vb(A)_1(t) vb(Phi)_1(t, sigma)$ and $vb(Phi)_1(t, t) = vb(I)$:
    $
        pdv(vb(Phi)_(1,2), t)
        & = vb(Phi)_1(t, t) vb(A)_(1,2)(t) vb(Phi)_2(t, tau) + integral_tau^t vb(A)_1(t) vb(Phi)_1(t, sigma) vb(A)_(1,2)(sigma) vb(Phi)_2(sigma, tau) dd(sigma) \
        & = vb(A)_(1,2)(t) vb(Phi)_2(t, tau) + vb(A)_1(t) vb(Phi)_(1,2)(t, tau)
    $
    as required. If $vb(A)_1$ is $n_1 times n_1$ and $vb(A)_2$ is $n_2 times n_2$, the integrand is $(n_1 times n_1)(n_1 times n_2)(n_2 times n_2) = n_1 times n_2$, so the formula is consistent even though $vb(A)_(1,2)$ need not be square. A numerical check on a sample system with a $2 times 1$ coupling block, comparing the integrated STM of the full system with the formula (Python):
    #raw(read("../outputs/text/s02.txt"), block: true)
<hwk:s02b>

#pagebreak(weak: true)

#prob[
    Given the state-space model:
    $
        dot(x)_1(t) & = t x_2(t) \
        dot(x)_2(t) & = -x_2(t) + e^(-t) u(t) \
               y(t) & = t x_1(t)
    $

    + Find the state transition matrix $Phi(t, tau)$.
    <hwk:p03a>

    + Find the impulse response $G(t, tau)$.
    <hwk:p03b>

    + Find the response $y(t)$ for $t gt.eq 1$ when $vb(x)(1)^TT = vecrow(1, -2)$ and $u(t) = 1$ for all $t gt.eq 1$.
    <hwk:p03c>
] <hwk:p03>

+ State transition matrix

    In matrix form, $dot(vb(x)) = A(t) vb(x) + B(t) u$ and $y = C(t) vb(x)$ with
    $
        A(t) = mat(0, t; 0, -1), quad B(t) = vec(0, e^(-t)), quad C(t) = mat(t, 0)
    $
    $A(t)$ is upper triangular, so the homogeneous equations can be solved from the bottom up, starting from $vb(x)(tau)$:
    $
        dot(x)_2 = -x_2 quad => quad x_2(t) = e^(-(t - tau)) x_2(tau)
    $
    $
        dot(x)_1 = t x_2(t) = t e^(-(t - tau)) x_2(tau)
        quad => quad
        x_1(t) = x_1(tau) + x_2(tau) integral_tau^t sigma e^(-(sigma - tau)) dd(sigma)
    $
    Using $integral sigma e^(-sigma) dd(sigma) = -(sigma + 1) e^(-sigma)$:
    $
        integral_tau^t sigma e^(-(sigma - tau)) dd(sigma)
        = e^tau lr([-(sigma + 1) e^(-sigma)])_tau^t
        = (tau + 1) - (t + 1) e^(-(t - tau))
    $
    The columns of $Phi(t, tau)$ are the responses to $vb(x)(tau) = vecrow(1, 0)^TT$ and $vb(x)(tau) = vecrow(0, 1)^TT$:
    $
        Phi(t, tau) = mat(1, (tau + 1) - (t + 1) e^(-(t - tau)); 0, e^(-(t - tau))) quad qed
    $
    Check: $Phi(tau, tau) = vb(I)$, and $pdv(Phi_(12), t) = t e^(-(t - tau)) = t Phi_(22)$, $pdv(Phi_(22), t) = -Phi_(22)$, as required by $pdv(Phi, t) = A(t) Phi$. This is also the scalar-block case of #link(<hwk:s02b>)[Problem 2] with $Phi_1 = 1$, $Phi_2 = e^(-(t - tau))$, and $A_(1,2) = t$. Compared with the numerically integrated STM (Python):
    #raw(read("../outputs/text/s03a.txt"), block: true)
<hwk:s03a>

+ Impulse response

    With $D = 0$, the impulse response is $G(t, tau) = C(t) Phi(t, tau) B(tau)$ for $t gt.eq tau$, and zero for $t < tau$:
    $
        C(t) Phi(t, tau) B(tau)
        = mat(t, 0) mat(1, (tau + 1) - (t + 1) e^(-(t - tau)); 0, e^(-(t - tau))) vec(0, e^(-tau))
        = t e^(-tau) [(tau + 1) - (t + 1) e^(-(t - tau))]
    $
    $
        G(t, tau) = cases(
            t [(tau + 1) e^(-tau) - (t + 1) e^(-t)] & quad t gt.eq tau,
            0 & quad t < tau,
        ) quad qed
    $
    An impulse at $tau$ sets $vb(x)(tau^+) = B(tau)$, so simulating the homogeneous system from that state reproduces $G(t, tau)$ (Python):
    #raw(read("../outputs/text/s03b.txt"), block: true)
<hwk:s03b>

+ Forced response from $t = 1$

    By superposition of the zero-input and zero-state responses,
    $
        y(t) = C(t) Phi(t, 1) vb(x)(1) + integral_1^t G(t, sigma) u(sigma) dd(sigma)
    $
    #emph[Zero-input part.] With $vb(x)(1) = vecrow(1, -2)^TT$:
    $
        Phi(t, 1) vb(x)(1)
        = vec(1 - 2[2 - (t + 1) e^(-(t - 1))], -2 e^(-(t - 1)))
        = vec(-3 + 2(t + 1) e^(1 - t), -2 e^(1 - t))
        quad => quad
        y_("zi")(t) = t [-3 + 2(t + 1) e^(1 - t)]
    $
    #emph[Zero-state part.] With $u = 1$ and $integral (sigma + 1) e^(-sigma) dd(sigma) = -(sigma + 2) e^(-sigma)$:
    $
        y_("zs")(t) & = t integral_1^t [(sigma + 1) e^(-sigma) - (t + 1) e^(-t)] dd(sigma)
                      = t [3 e^(-1) - (t + 2) e^(-t) - (t + 1)(t - 1) e^(-t)] \
                    & = t [3 e^(-1) - (t^2 + t + 1) e^(-t)]
    $
    Adding the two parts and collecting the $e^(-t)$ terms:
    $
        y(t) = t [3/e - 3 + (2e(t + 1) - t^2 - t - 1) e^(-t)], quad t gt.eq 1 quad qed
    $
    Check: $y(1) = 3/e - 3 + (4e - 3)/e = 1 = C(1) vb(x)(1)$. As $t -> oo$, $x_1(t) -> 3/e - 3$, so the output grows without bound like $y(t) approx (3/e - 3) t approx -1.896 t$ because of the time-varying output gain $C(t) = mat(t, 0)$. The closed form matches a direct simulation of the state equations (Python):
    #raw(read("../outputs/text/s03c.txt"), block: true)
    #figure(
        image("../outputs/figures/s03c.png", width: 90%),
        caption: [Response $y(t)$ for $vb(x)(1)^TT = vecrow(1, -2)$ and $u = 1$, with its zero-input and zero-state parts.],
    )
<hwk:s03c>

#pagebreak(weak: true)

#prob[
    The characteristic polynomial of the $vb(A)$ matrix for an LTI dynamic system can be written as $r(lambda) = (lambda^2 + 2 lambda + 10)(lambda + 3)^3$. You also know that the rank of $(vb(A) + 3 vb(I))$ is 4.

    Let $vb(P)$ be the matrix which will transform $vb(A)$ to an appropriate real block modal form. Assume the $vb(P)$ matrix is such that the components of the modal form arising from the complex eigenvalues will be located in the upper left corner.

    + Compute the matrix $vb(P)^(-1) e^(vb(A) t) vb(P)$.
    <hwk:p04a>

    + Suppose that we identify the outputs of this system to be $vb(y)(t) = vb(P)^(-1) vbu(x)(t)$. If the initial condition is $vbu(x)_0 = vbu(p)_4$ where $vbu(p)_4$ is the fourth column of the $vb(P)$ matrix, find $vbu(y)(t)$ explicitly as a function of time.
    <hwk:p04b>

    + Suppose in #link(<hwk:p04b>)[Part B] that the initial condition is instead $vbu(x)_0 = vbu(p)_2$. Find $vbu(y)(t)$.
    <hwk:p04c>
] <hwk:p04>

+ Modal form and $vb(P)^(-1) e^(vb(A) t) vb(P)$

    $r(lambda)$ has degree 5, so $vb(A)$ is $5 times 5$. Its eigenvalues are $lambda = -1 plus.minus 3j$ (from $lambda^2 + 2 lambda + 10 = (lambda + 1)^2 + 9$) and $lambda = -3$ with algebraic multiplicity 3. The geometric multiplicity of $lambda = -3$ is
    $
        dim cal(N)(vb(A) + 3 vb(I)) = 5 - rank(vb(A) + 3 vb(I)) = 5 - 4 = 1
    $
    so $lambda = -3$ has only one independent eigenvector, and its three-fold eigenvalue forms a single $3 times 3$ Jordan block. (Two or three blocks would give $rank(vb(A) + 3 vb(I)) = 3$ or $2$.) The columns of $vb(P)$ are therefore:
    - $vbu(p)_1 = Re(vbu(v))$ and $vbu(p)_2 = Im(vbu(v))$, where $vbu(v)$ is an eigenvector for $-1 + 3j$. As in #link(<hwk:s01a>)[Problem 1],
        $
            vb(A) mat(vbu(p)_1, vbu(p)_2) = mat(vbu(p)_1, vbu(p)_2) mat(-1, 3; -3, -1)
        $
    - A chain of generalized eigenvectors for $lambda = -3$:
        $
            (vb(A) + 3 vb(I)) vbu(p)_3 = vb(0), quad
            (vb(A) + 3 vb(I)) vbu(p)_4 = vbu(p)_3, quad
            (vb(A) + 3 vb(I)) vbu(p)_5 = vbu(p)_4
        $
        Equivalently, $vb(A) vbu(p)_3 = -3 vbu(p)_3$, $vb(A) vbu(p)_4 = vbu(p)_3 - 3 vbu(p)_4$, and $vb(A) vbu(p)_5 = vbu(p)_4 - 3 vbu(p)_5$.

    Together these give the real block-modal (real Jordan) form:
    $
        vb(J) = vb(P)^(-1) vb(A) vb(P) = mat(
            -1, 3, 0, 0, 0;
            -3, -1, 0, 0, 0;
            0, 0, -3, 1, 0;
            0, 0, 0, -3, 1;
            0, 0, 0, 0, -3;
        )
    $
    Since $vb(P)^(-1) vb(A)^k vb(P) = vb(J)^k$, the exponential series gives $vb(P)^(-1) e^(vb(A) t) vb(P) = e^(vb(J) t)$, which can be computed block by block:
    - The complex block is $-vb(I) + 3 vb(S)$ with $vb(S) = mat(0, 1; -1, 0)$. These two terms commute, and $vb(S)^2 = -vb(I)$ gives $e^(theta vb(S)) = cos(theta) vb(I) + sin(theta) vb(S)$, so
        $
            e^((-vb(I) + 3 vb(S)) t) = e^(-t) e^(3t vb(S)) = e^(-t) mat(cos(3t), sin(3t); -sin(3t), cos(3t))
        $
    - The Jordan block is $-3 vb(I) + vb(N)$, where $vb(N)$ has ones on the superdiagonal, so $vb(N)^3 = vb(0)$. These two terms also commute, so the series for $e^(vb(N) t)$ terminates:
        $
            e^((-3 vb(I) + vb(N)) t) = e^(-3t) (vb(I) + t vb(N) + t^2/2 vb(N)^2) = e^(-3t) mat(1, t, t^2/2; 0, 1, t; 0, 0, 1)
        $

    Thus
    $
        vb(P)^(-1) e^(vb(A) t) vb(P) = e^(vb(J) t) = mat(
            e^(-t) cos(3t), e^(-t) sin(3t), 0, 0, 0;
            -e^(-t) sin(3t), e^(-t) cos(3t), 0, 0, 0;
            0, 0, e^(-3t), t e^(-3t), t^2/2 e^(-3t);
            0, 0, 0, e^(-3t), t e^(-3t);
            0, 0, 0, 0, e^(-3t);
        ) quad qed
    $
<hwk:s04a>

+ Response from $vbu(x)_0 = vbu(p)_4$

    With $vbu(x)(t) = e^(vb(A) t) vbu(x)_0$:
    $
        vbu(y)(t) = vb(P)^(-1) e^(vb(A) t) vbu(x)_0 = (vb(P)^(-1) e^(vb(A) t) vb(P)) vb(P)^(-1) vbu(x)_0 = e^(vb(J) t) vb(P)^(-1) vbu(x)_0
    $
    Since $vbu(p)_4$ is the fourth column of $vb(P)$, $vb(P)^(-1) vbu(p)_4 = vbu(e)_4$, so $vbu(y)(t)$ is the fourth column of $e^(vb(J) t)$:
    $
        vbu(y)(t) = vec(0, 0, t e^(-3t), e^(-3t), 0) quad qed
    $
    In the original coordinates, $vbu(x)(t) = e^(-3t) (t vbu(p)_3 + vbu(p)_4)$. Starting on the generalized eigenvector $vbu(p)_4$ also drives the eigenvector direction $vbu(p)_3$ through the secular term $t e^(-3t)$. It never excites $vbu(p)_5$ or the oscillatory mode.
<hwk:s04b>

+ Response from $vbu(x)_0 = vbu(p)_2$

    Now $vb(P)^(-1) vbu(p)_2 = vbu(e)_2$, so $vbu(y)(t)$ is the second column of $e^(vb(J) t)$:
    $
        vbu(y)(t) = vec(e^(-t) sin(3t), e^(-t) cos(3t), 0, 0, 0) quad qed
    $
    In the original coordinates, $vbu(x)(t) = e^(-t) (sin(3t) vbu(p)_1 + cos(3t) vbu(p)_2)$, a decaying oscillation confined to the oscillatory mode.

    As a numerical check, $vb(A) = vb(P) vb(J) vb(P)^(-1)$ was built from a random invertible $vb(P)$. The Python calculation confirms that its characteristic polynomial and ranks match the given data, and that all three results above match $vb(P)^(-1) e^(vb(A) t)$ computed with a matrix exponential:
    #raw(read("../outputs/text/s04.txt"), block: true)
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
            vb(Phi)(t, tau) = e^(vb(A)_1 integral_tau^t f_1(sigma) dd(sigma)) e^(vb(A)_2 integral_tau^t f_2(sigma) dd(sigma)) space ... space e^(vb(A)_r integral_tau^t f_r(sigma) dd(sigma))
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

+ STM as a matrix exponential

    Define
    $
        F_j (t, tau) = integral_tau^t f_j (sigma) dd(sigma), quad
        vb(M)(t, tau) = integral_tau^t vb(A)(sigma) dd(sigma) = sum_(j=1)^r F_j (t, tau) vb(A)_j
    $
    Because the $vb(A)_j$ commute pairwise, $vb(A)(t)$ commutes with $vb(M)(t, tau)$:
    $
        vb(A)(t) vb(M)(t, tau)
        = sum_(i=1)^r sum_(j=1)^r f_i (t) F_j (t, tau) vb(A)_i vb(A)_j
        = sum_(i=1)^r sum_(j=1)^r F_j (t, tau) f_i (t) vb(A)_j vb(A)_i
        = vb(M)(t, tau) vb(A)(t)
    $
    Since $pdv(vb(M), t) = vb(A)(t)$ commutes with $vb(M)$, the product rule collapses:
    $
        pdv(, t) vb(M)^k = sum_(m=0)^(k-1) vb(M)^m vb(A)(t) vb(M)^(k-1-m) = k vb(A)(t) vb(M)^(k-1)
    $
    Assume the $f_j$ are continuous, so $vb(A)(t)$ and $vb(M)(t, tau)$ are bounded on any bounded $t$-interval. Then the differentiated series $vb(A)(t) sum_(k >= 1) vb(M)^(k-1)/(k-1)!$ converges uniformly there, which justifies differentiating the exponential series term by term:
    $
        pdv(, t) e^(vb(M)(t, tau))
        = sum_(k=1)^oo k/k! vb(A)(t) vb(M)^(k-1)
        = vb(A)(t) sum_(k=1)^oo vb(M)^(k-1)/(k-1)!
        = vb(A)(t) e^(vb(M)(t, tau))
    $
    Also $e^(vb(M)(tau, tau)) = e^(vb(0)) = vb(I)$. So $e^(vb(M)(t, tau))$ solves the defining problem for the STM, whose solution is unique:
    $
        vb(Phi)(t, tau) = e^(integral_tau^t vb(A)(sigma) dd(sigma)) = exp(sum_(j=1)^r vb(A)_j integral_tau^t f_j (sigma) dd(sigma)) quad qed
    $
    Without commutativity, the derivative of $vb(M)^k$ does not collapse to $k vb(A) vb(M)^(k-1)$, and this formula generally fails (see #link(<hwk:s06c>)[Problem 6, Part C]).
<hwk:s05a>

+ Factoring the exponential

    Let $vb(X)_j = F_j (t, tau) vb(A)_j$. These commute pairwise, since $vb(X)_i vb(X)_j = F_i F_j vb(A)_i vb(A)_j = F_j F_i vb(A)_j vb(A)_i = vb(X)_j vb(X)_i$.

    For commuting $vb(X)$ and $vb(Y)$, the binomial theorem holds: $(vb(X) + vb(Y))^k = sum_(m=0)^k binom(k, m) vb(X)^m vb(Y)^(k-m)$. Multiplying the two absolutely convergent series as a Cauchy product,
    $
        e^(vb(X) + vb(Y))
        = sum_(k=0)^oo sum_(m=0)^k (vb(X)^m vb(Y)^(k-m))/(m! (k-m)!)
        = (sum_(m=0)^oo vb(X)^m/m!) (sum_(n=0)^oo vb(Y)^n/n!)
        = e^(vb(X)) e^(vb(Y))
    $
    $vb(X)_1$ commutes with $vb(X)_2 + dots.c + vb(X)_r$, so applying this repeatedly (induction on $r$):
    $
        e^(vb(X)_1 + vb(X)_2 + dots.c + vb(X)_r)
        = e^(vb(X)_1) e^(vb(X)_2 + dots.c + vb(X)_r)
        = dots.c
        = e^(vb(X)_1) e^(vb(X)_2) dots.c e^(vb(X)_r)
    $
    Combining with #link(<hwk:s05a>)[Part A]:
    $
        vb(Phi)(t, tau) = e^(vb(A)_1 integral_tau^t f_1(sigma) dd(sigma)) e^(vb(A)_2 integral_tau^t f_2(sigma) dd(sigma)) dots.c e^(vb(A)_r integral_tau^t f_r(sigma) dd(sigma)) quad qed
    $
    The factors commute, so their order does not matter. As a numerical check, three commuting matrices (polynomials in one random matrix) were used with $f_j = cos(t), t, e^(-t)$, and both forms were compared with the numerically integrated STM (Python):
    #raw(read("../outputs/text/s05b.txt"), block: true)
<hwk:s05b>

+ STM for the rotation-like $vb(A)(t)$

    Split $vb(A)(t)$ into constant matrices:
    $
        vb(A)(t) = cos(3t) vb(I) + sin(3t) vb(S), quad vb(S) = mat(0, 1; -1, 0)
    $
    This has the required form with $r = 2$: $f_1 = cos(3t)$, $vb(A)_1 = vb(I)$, $f_2 = sin(3t)$, and $vb(A)_2 = vb(S)$. Since $vb(I)$ commutes with every matrix, the results above apply, with
    $
        F_1 = integral_tau^t cos(3 sigma) dd(sigma) = (sin(3t) - sin(3 tau))/3, quad
        F_2 = integral_tau^t sin(3 sigma) dd(sigma) = (cos(3 tau) - cos(3t))/3
    $
    By #link(<hwk:s05b>)[Part B], $vb(Phi) = e^(F_1 vb(I)) e^(F_2 vb(S))$. First, $e^(F_1 vb(I)) = e^(F_1) vb(I)$. Second, since $vb(S)^2 = -vb(I)$, the series for $e^(theta vb(S))$ splits into even powers $(-1)^k vb(I)$ and odd powers $(-1)^k vb(S)$:
    $
        e^(theta vb(S)) = cos(theta) vb(I) + sin(theta) vb(S) = mat(cos(theta), sin(theta); -sin(theta), cos(theta))
    $
    Therefore
    $
        vb(Phi)(t, tau) = e^((sin(3t) - sin(3 tau))\/3) mat(cos(theta), sin(theta); -sin(theta), cos(theta)), quad
        theta(t, tau) = (cos(3 tau) - cos(3t))/3 quad qed
    $
    Check: $vb(Phi)(tau, tau) = vb(I)$, and $pdv(vb(Phi), t) = cos(3t) vb(Phi) + sin(3t) vb(S) vb(Phi) = vb(A)(t) vb(Phi)$. Python confirms the commutativity and compares the result with the numerically integrated STM:
    #raw(read("../outputs/text/s05c.txt"), block: true)
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

+ Diagonal $vb(A)(t)$

    The equations decouple into two scalar equations, each solved by separation of variables:
    $
        dot(x)_1 = t x_1 quad => quad x_1(t) = e^(integral_tau^t sigma dd(sigma)) x_1(tau) = e^((t^2 - tau^2)\/2) x_1(tau)
    $
    $
        dot(x)_2 = cos(t) x_2 quad => quad x_2(t) = e^(integral_tau^t cos(sigma) dd(sigma)) x_2(tau) = e^(sin(t) - sin(tau)) x_2(tau)
    $
    $
        vb(Phi)(t, tau) = mat(e^((t^2 - tau^2)\/2), 0; 0, e^(sin(t) - sin(tau))) quad qed
    $
    This agrees with #link(<hwk:s05a>)[Problem 5]: diagonal matrices commute, so $vb(Phi) = e^(integral_tau^t vb(A)(sigma) dd(sigma))$.
<hwk:s06a>

+ Rank-one $vb(A)(t)$

    The second equation is $dot(x)_2 = 0$, so $x_2(t) = x_2(tau)$. The first equation is $dot(x)_1 = cos(2t) (x_1 + x_2)$. With $w = x_1 + x_2$ and $dot(w) = dot(x)_1$:
    $
        dot(w) = cos(2t) w
        quad => quad
        w(t) = e^(F(t, tau)) w(tau), quad
        F(t, tau) = integral_tau^t cos(2 sigma) dd(sigma) = (sin(2t) - sin(2 tau))/2
    $
    $
        x_1(t) = w(t) - x_2(t) = e^F x_1(tau) + (e^F - 1) x_2(tau)
    $
    $
        vb(Phi)(t, tau) = mat(e^(F(t, tau)), e^(F(t, tau)) - 1; 0, 1), quad F(t, tau) = (sin(2t) - sin(2 tau))/2 quad qed
    $
    Cross-check with #link(<hwk:s05a>)[Problem 5]: $vb(A)(t) = cos(2t) vb(N)$ uses a single constant matrix (so it trivially commutes), and $vb(N)$ is idempotent:
    $
        vb(N) = mat(1, 1; 0, 0), quad vb(N)^2 = vb(N)
        quad => quad
        e^(F vb(N)) = vb(I) + sum_(k=1)^oo F^k/k! vb(N) = vb(I) + (e^F - 1) vb(N)
    $
    which is the same matrix.
<hwk:s06b>

+ Upper-triangular $vb(A)(t)$ with time-varying coupling

    Here $vb(A)(t)$ does not commute with itself at different times:
    $
        vb(A)(t) vb(A)(s) = mat(4, -2 sin(s) - sin(t); 0, 1)
        != mat(4, -2 sin(t) - sin(s); 0, 1) = vb(A)(s) vb(A)(t)
        quad "unless" sin(t) = sin(s)
    $
    so the exponential of $integral vb(A)$ is #emph[not] the STM. Instead, $vb(A)(t)$ is upper triangular, so #link(<hwk:s02b>)[Problem 2] applies with scalar blocks $A_1 = -2$, $A_2 = -1$, and $A_(1,2) = sin(t)$:
    $
        Phi_1(t, tau) = e^(-2(t - tau)), quad
        Phi_2(t, tau) = e^(-(t - tau))
    $
    $
        Phi_(1,2)(t, tau)
        = integral_tau^t e^(-2(t - sigma)) sin(sigma) e^(-(sigma - tau)) dd(sigma)
        = e^(-2t + tau) integral_tau^t e^sigma sin(sigma) dd(sigma)
    $
    Using $integral e^sigma sin(sigma) dd(sigma) = 1/2 e^sigma (sin(sigma) - cos(sigma))$:
    $
        Phi_(1,2)(t, tau)
        = 1/2 e^(-2t + tau) [e^t (sin(t) - cos(t)) - e^tau (sin(tau) - cos(tau))]
        = 1/2 [e^(-(t - tau)) (sin(t) - cos(t)) - e^(-2(t - tau)) (sin(tau) - cos(tau))]
    $
    $
        vb(Phi)(t, tau) = mat(
            e^(-2(t - tau)), 1/2 [e^(-(t - tau)) (sin(t) - cos(t)) - e^(-2(t - tau)) (sin(tau) - cos(tau))];
            0, e^(-(t - tau));
        ) quad qed
    $
    Check: $Phi_(1,2)(tau, tau) = 0$, and
    $
        pdv(Phi_(1,2), t) = e^(-(t - tau)) cos(t) + e^(-2(t - tau)) (sin(tau) - cos(tau)) = -2 Phi_(1,2) + sin(t) e^(-(t - tau))
    $
    which is the $(1, 2)$ entry of $vb(A)(t) vb(Phi)$.
<hwk:s06c>

+ Dyadic expansion for constant $vb(A)$

    For a constant, diagonalizable $vb(A) = vb(V) vb(Lambda) vb(V)^(-1)$, let $vbu(v)_i$ be the eigenvectors (columns of $vb(V)$) and $vbu(w)_i^TT$ the reciprocal basis (rows of $vb(V)^(-1)$, so $vbu(w)_i^TT vbu(v)_j = delta_(i j)$). Then $vb(A)$ has the dyadic expansion
    $
        vb(A) = sum_i lambda_i vbu(v)_i vbu(w)_i^TT
        quad => quad
        vb(A)^k = sum_i lambda_i^k vbu(v)_i vbu(w)_i^TT
        quad => quad
        vb(Phi)(t, tau) = e^(vb(A)(t - tau)) = sum_i e^(lambda_i (t - tau)) vbu(v)_i vbu(w)_i^TT
    $
    The eigenvalues satisfy $det(lambda vb(I) - vb(A)) = (lambda + 3)^2 - 1 = (lambda + 2)(lambda + 4) = 0$, so $lambda_1 = -2$ and $lambda_2 = -4$. Eigenvectors:
    $
        (vb(A) + 2 vb(I)) vbu(v)_1 = mat(-1, -1; -1, -1) vbu(v)_1 = vb(0) => vbu(v)_1 = vec(1, -1), quad
        (vb(A) + 4 vb(I)) vbu(v)_2 = mat(1, -1; -1, 1) vbu(v)_2 = vb(0) => vbu(v)_2 = vec(1, 1)
    $
    Reciprocal basis:
    $
        vb(V)^(-1) = mat(1, 1; -1, 1)^(-1) = 1/2 mat(1, -1; 1, 1)
        quad => quad
        vbu(w)_1^TT = 1/2 vecrow(1, -1), quad
        vbu(w)_2^TT = 1/2 vecrow(1, 1)
    $
    The dyads are
    $
        vbu(v)_1 vbu(w)_1^TT = 1/2 mat(1, -1; -1, 1), quad
        vbu(v)_2 vbu(w)_2^TT = 1/2 mat(1, 1; 1, 1)
    $
    They sum to $vb(I)$, and $-2 vbu(v)_1 vbu(w)_1^TT - 4 vbu(v)_2 vbu(w)_2^TT = vb(A)$, as expected. Therefore
    $
        vb(Phi)(t, tau)
        = e^(-2(t - tau)) / 2 mat(1, -1; -1, 1) + e^(-4(t - tau)) / 2 mat(1, 1; 1, 1)
        = 1/2 mat(
            e^(-2(t - tau)) + e^(-4(t - tau)), -e^(-2(t - tau)) + e^(-4(t - tau));
            -e^(-2(t - tau)) + e^(-4(t - tau)), e^(-2(t - tau)) + e^(-4(t - tau));
        ) quad qed
    $
    Each STM above was compared with the numerically integrated STM. For Part C, Python also shows that $exp(integral vb(A))$ is not the STM, and for Part D it compares the dyadic sum with a matrix exponential:
    #raw(read("../outputs/text/s06.txt"), block: true)
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

+ Frozen-time eigenvalues

    Abbreviate $s = sin(t)$ and $c = cos(t)$. The trace and determinant of $vb(A)(t)$ are
    $
        trace vb(A)(t) = -2 + 1.5(c^2 + s^2) = -1/2
    $
    $
        det vb(A)(t) & = (-1 + 1.5 c^2)(-1 + 1.5 s^2) - (1 - 1.5 s c)(-1 - 1.5 s c) \
                     & = (1 - 1.5(c^2 + s^2) + 2.25 s^2 c^2) + (1 - 2.25 s^2 c^2) = 1/2
    $
    Both are independent of $t$, so the characteristic polynomial $lambda^2 - trace(vb(A)) lambda + det(vb(A)) = lambda^2 + 1/2 lambda + 1/2$ is the same at every instant:
    $
        lambda_(1,2) = (-1/2 plus.minus sqrt(1/4 - 2))/2 = -1/4 plus.minus j sqrt(7)/4 quad qed
    $
    The eigenvalues do not change with time, and $Re(lambda_(1,2)) = -1/4 < 0$, so they always lie in the open left half-plane. Python confirms this on a dense grid over a period:
    #raw(read("../outputs/text/s07a.txt"), block: true)
<hwk:s07a>

+ Fundamental solutions

    Differentiating $vb(psi)_1$ and computing $vb(A) vb(psi)_1$ (again with $s = sin(t)$, $c = cos(t)$):
    $
        dot(vb(psi))_1 = e^(-t) vec(c - s, -s - c), quad
        vb(A) vb(psi)_1
        = e^(-t) vec((-1 + 1.5 c^2) s + (1 - 1.5 s c) c, (-1 - 1.5 s c) s + (-1 + 1.5 s^2) c)
        = e^(-t) vec(c - s, -s - c)
    $
    Similarly for $vb(psi)_2$, using $c^2 + s^2 = 1$:
    $
        dot(vb(psi))_2 = e^(t\/2) vec(1/2 c - s, -1/2 s - c), quad
        vb(A) vb(psi)_2
        = e^(t\/2) vec((-1 + 1.5 c^2) c - (1 - 1.5 s c) s, (-1 - 1.5 s c) c - (-1 + 1.5 s^2) s)
        = e^(t\/2) vec(-c + 1.5 c - s, -c + s - 1.5 s)
    $
    which equals $dot(vb(psi))_2$. Both are therefore solutions of $dot(vbu(x)) = vb(A)(t) vbu(x)$. They are linearly independent for every $t$:
    $
        det mat(vb(psi)_1, vb(psi)_2) = e^(-t) e^(t\/2) det mat(s, c; c, -s) = -e^(-t\/2) != 0 quad qed
    $
    This is consistent with the Abel--Jacobi--Liouville formula, $det vb(Psi)(t) = det vb(Psi)(0) e^(integral_0^t trace vb(A)) = -e^(-t\/2)$. Hence $vb(psi)_1$ and $vb(psi)_2$ form a fundamental set. Python confirms both properties:
    #raw(read("../outputs/text/s07b.txt"), block: true)
<hwk:s07b>

+ State transition matrix and stability

    Collect the fundamental solutions into the fundamental matrix $vb(Psi)(t)$, whose columns are $vb(psi)_1(t)$ and $vb(psi)_2(t)$. With $det vb(Psi)(tau) = -e^(-tau\/2)$:
    $
        vb(Psi)(t) = mat(e^(-t) sin(t), e^(t\/2) cos(t); e^(-t) cos(t), -e^(t\/2) sin(t)), quad
        vb(Psi)^(-1)(tau) = mat(e^tau sin(tau), e^tau cos(tau); e^(-tau\/2) cos(tau), -e^(-tau\/2) sin(tau))
    $
    Then $vb(Phi)(t, tau) = vb(Psi)(t) vb(Psi)^(-1)(tau)$:
    $
        vb(Phi)(t, tau) = mat(
            e^(-(t - tau)) sin(t) sin(tau) + e^((t - tau)\/2) cos(t) cos(tau), e^(-(t - tau)) sin(t) cos(tau) - e^((t - tau)\/2) cos(t) sin(tau);
            e^(-(t - tau)) cos(t) sin(tau) - e^((t - tau)\/2) sin(t) cos(tau), e^(-(t - tau)) cos(t) cos(tau) + e^((t - tau)\/2) sin(t) sin(tau);
        ) quad qed
    $
    Clearly $vb(Phi)(tau, tau) = vb(I)$. Multiplying out confirms the equivalent factored form
    $
        vb(Phi)(t, tau) = vb(R)(t) mat(e^((t - tau)\/2), 0; 0, e^(-(t - tau))) vb(R)(tau)^TT, quad
        vb(R)(t) = mat(cos(t), sin(t); -sin(t), cos(t))
    $
    where $vb(R)(t)$ is a rotation matrix.

    #emph[Stability.] The system is #strong[not] stable. Starting from $vbu(x)_0 = vecrow(1, 0)^TT$ at $t_0 = 0$ gives
    $
        vbu(x)(t) = vb(Phi)(t, 0) vec(1, 0) = vb(psi)_2(t) = e^(t\/2) vec(cos(t), -sin(t)),
        quad
        norm(vbu(x)(t)) = e^(t\/2) -> oo quad qed
    $
    In general, $vbu(x)(t) = vb(Psi)(t) vb(Psi)^(-1)(t_0) vbu(x)_0$ decays only if the $vb(psi)_2$ coefficient $e^(-t_0\/2)(cos(t_0) x_(0,1) - sin(t_0) x_(0,2))$ is zero. That is, $vbu(x)_0$ must be parallel to $vb(psi)_1(t_0)$, a one-dimensional set. Every other initial state diverges like $e^(t\/2)$.

    The same conclusion follows from Floquet theory: $vb(A)(t)$ has period $pi$, and the monodromy matrix is $vb(Phi)(pi, 0) = op("diag")(-e^(pi\/2), -e^(-pi))$. One Floquet multiplier has magnitude $e^(pi\/2) approx 4.81 > 1$.
    #raw(read("../outputs/text/s07c.txt"), block: true)
    #figure(
        image("../outputs/figures/s07c.png", width: 90%),
        caption: [State norm for three initial states, compared with the $e^(-t\/4)$ decay suggested by the frozen-time eigenvalues.],
    )
<hwk:s07c>

+ Eigenvalues of a time-varying $vb(A)(t)$

    In #link(<hwk:s07a>)[Part A], the eigenvalues of $vb(A)(t)$ are constant and in the left half-plane for every $t$. Read as if the system were LTI, they would predict decay like $e^(-t\/4)$. Yet #link(<hwk:s07c>)[Part C] shows that almost every solution grows like $e^(t\/2)$. For a time-varying system, the instantaneous ("frozen-time") eigenvalues of $vb(A)(t)$ do not determine stability:
    - Left-half-plane frozen eigenvalues are not sufficient for stability, as this example shows. Neither are they necessary: frozen eigenvalues in the right half-plane can coexist with asymptotic stability (an example follows below).
    - The reason is visible here. $vb(A)(t) = vb(R)(t) vb(A)(0) vb(R)(t)^TT$ is a rotating copy of $vb(A)(0)$, with $vb(R)(t)$ from #link(<hwk:s07c>)[Part C]. A similarity transformation preserves the frozen eigenvalues, but it ignores the motion of the coordinate frame. In the rotating coordinates $vbu(z) = vb(R)(t)^TT vbu(x)$, using $dot(vb(R)) = vb(S) vb(R)$, the dynamics become time-invariant:
        $
            dot(vbu(z)) = (vb(A)(0) - vb(S)) vbu(z), quad
            vb(A)(0) = mat(1/2, 1; -1, -1), quad
            vb(S) = mat(0, 1; -1, 0), quad
            vb(A)(0) - vb(S) = mat(1/2, 0; 0, -1)
        $
        This has an unstable eigenvalue $+1/2$. The frame rotation contributes $-vb(S)$, which the frozen eigenvalues never see. The same mechanism can also stabilize a system. Rotating $vb(A)_0 = mat(-1, 1; 9, -1)$ instead (frozen eigenvalues $2$ and $-4$) gives $vb(A)_0 - vb(S) = mat(-1, 0; 10, -1)$, which is Hurwitz. That system is asymptotically stable even though one frozen eigenvalue is in the right half-plane for every $t$.
    - Stability of an LTV system must be judged from $vb(Phi)(t, tau)$ itself (boundedness and decay), from Lyapunov arguments, or, for periodic systems, from the Floquet multipliers (the eigenvalues of the monodromy matrix $vb(Phi)(T, 0)$). Frozen-time eigenvalues are only trustworthy under extra conditions. One case is when $vb(A)(t)$ is bounded, commutes with itself at all times (as in #link(<hwk:s05a>)[Problem 5]), #emph[and] its eigenvalues satisfy $Re(lambda_i (t)) <= -epsilon < 0$ uniformly. Another is when $vb(A)(t)$ is bounded, satisfies the same uniform bound $Re(lambda_i (t)) <= -epsilon$, and varies slowly enough compared with that decay rate. Commutativity alone is not enough: $dot(x) = -x\/(1 + t)^2$ has a negative frozen eigenvalue for all $t$, yet $x(t) = x_0 e^(1\/(1 + t) - 1) -> x_0\/e != 0$.

    So the eigenvalues of a time-varying $vb(A)(t)$ are not a reliable stability test. $quad qed$
<hwk:s07d>

#pagebreak(weak: true)

== Code

#codly(header: [./src/index.py])
#raw(read("../src/index.py"), block: true, lang: "python") <code:index.py>
