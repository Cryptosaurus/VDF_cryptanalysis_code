
This code is an alternative implementation of the complexity estimates. It uses the following approximations.

## Additional formulas 

### Almost smoothness

Let $\Psi(x,y)$ be the number of $y$-smooth numbers less or equal to $x$.

The smoothness probability $P_{sm,x}(y)$ can then be defined as $\frac{\Psi(x,y)}{x}$ and is known to be approximated.

Let $\Psi^A(x,y,z)$ be the number of numbers less or equal to $x$ such that, apart from a single factor $q<z$, all factors are smaller than $y$.

We have effectively that every almost-smooth number is product of a  prime $p$ between $y$ and $z$ and a smooth number below $x/p$.
Thus we have
$$
\Psi^A(x,y,z) = \sum_{\text{prime }p\in[y;z)}\Psi(x/p,y).
$$
And the almost-smoothness probability is defined as 
$$
P_{sm,x}(y,z) = \frac{\Psi^A(x,y,z)}{x} =  \sum_{\text{prime }p\in[y;z)} \frac{\Psi(x/p,y)}{x} =  \sum_{\text{prime }p\in[y;z)} \frac{P_{sm,x/p}(y)}{p}
$$
For $y,z<<x$ we assume that  $P_{sm,x/p}(y)$ is the same for all $p$ of the same bitlength $l$. Switching to powers of two for $x,y,z$ we have
\begin{align}
P_{sm,2^x}^A(2^y,2^z) = &\sum_{\text{prime }p\in[2^y;2^z)} \frac{P_{sm,2^{x-\log p}}(2^y)}{p}\approx
\sum_{i\in[y;z)} \frac{P_{sm,2^{x-i}}(2^y)}{2^i} \times
\{\text{prime }p\in[2^i;2^{i+1}]\}= \\
=&\sum_{i\in[y;z)} \frac{P_{sm,2^{x-i}}(2^y)}{2^i}(\frac{2^{i+1}}{(i+1)\ln 2}-\frac{2^{i}}{i\ln 2}) = 
\sum_{i\in[y;z)}\frac{P_{sm,2^{x-i}}(2^y)(i-1)}{i(i+1)\ln 2}
\end{align}
Given an approximation for $P_{sm,x}$ we can approximate $P_{sm,2^x}^A(2^y,2^z)$ via the sum above.

### Size of a smooth factor.

Consider the probability $P_{psm,2^x}(2^y,2^z)$ that the random number below $2^x$ has its $2^y$-smoothness part bigger than $2^z$. 

Note that all such numbers are product of smooth numbers and numbers without a factor in the smoothness area. Thus we have
$$
\Psi^P(2^x,2^y,2^z) = \sum_{a<2^{x-z} \text{ has no factors  }<2^y}\Psi(2^{x}/a,2^y).
$$
By Mertens' third theorem, the probability to have no factor below $2^y$ is approximated as $\frac{1}{e^{\gamma}y\ln 2}$. Thus we have the following approximation
\begin{align}
P_{psm,2^x}(2^y,2^z) =& \Psi^P(2^x,2^y,2^z)/2^x = 
\frac{1}{e^{\gamma}y\ln 2}\sum_{a<2^{x-z} }\frac{\Psi(2^{x}/a,2^y)}{2^x}= \\=&\frac{1}{e^{\gamma}y\ln 2}\sum_{a<2^{x-z} }\frac{P_{sm,2^{x-\log a}}(2^y))}{a}\approx\frac{1}{e^{\gamma}y\ln 2}\sum_{0\leq i<{x-z} }\frac{P_{sm,2^{x-i}}(2^y))}{\ln(2)}
\\
\approx&\frac{2^{-0.23}}{y}\sum_{i<{x-z} }P_{sm,2^{x-i}}(2^y)
\end{align}


### Sum of prime powers inversed

Well known that
$$
\sum_{p\leq n}\frac{1}{p}\approx \ln\ln n + 0.26
$$
and
$$
\sum_{\text{prime }p}\frac{1}{p^k} = P(k)
$$
where $P$ is the prime zeta function. Given that $P(k+1)<P(k)/2$ we get that
$$
\sum_{\text{prime }p, \text{ integer }k\geq 2}\frac{1}{p^k} \leq 2P(2) < 0.92
$$
and 
$$
\sum_{\text{prime }p, \text{ integer }k\geq 1}\frac{1}{p^k} \leq\ln\ln n + 1.2
$$