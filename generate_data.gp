\\ =========================================================================
\\ Memory and precision configuration
\\ =========================================================================
default(parisize, "2G");
default(realprecision, 30);   \\ reduced to speed up

\\ =========================================================================
\\ Definition of the L-functions
\\ =========================================================================
L = lfuncreate(Mod(3,4));    \\ L(s, chi_-4)
Z = lfuncreate(1);            \\ Riemann zeta

\\ =========================================================================
\\ Maximum height
\\ =========================================================================
Tmax = 10050;

\\ -------------------------------------------------------------------------
\\ Computation of the zeros
\\ -------------------------------------------------------------------------
print("Computing zeros of L(s, chi_-4) up to height ", Tmax, "...");
zerosL = lfunzeros(L, [0, Tmax]);

print("Computing zeros of Riemann zeta up to height ", Tmax, "...");
zerosZeta = lfunzeros(Z, [0, Tmax]);

print("Number of L(s,chi_-4) zeros: ", #zerosL);
print("Number of zeta zeros:       ", #zerosZeta);
print("First L zero: ", zerosL[1]);
print("Last  L zero: ", zerosL[#zerosL]);
print("First zeta zero: ", zerosZeta[1]);
print("Last  zeta zero: ", zerosZeta[#zerosZeta]);

\\ Save the raw zeros for potential reuse
write("zeros_L4_raw.txt", Str(zerosL));
write("zeros_zeta_raw.txt", Str(zerosZeta));
print("Raw zeros saved.");

\\ =========================================================================
\\ BLOCK 1: Zeros of L(s, chi_-4)
\\ For these zeros:  zeta_Q(i)'(rho) = zeta(rho) * L'(rho)
\\ =========================================================================
linha1 = List();
linha2 = List();

print("Processing L(s,chi_-4) zeros...");
for(i=1, #zerosL, {
    if(i % 500 == 0, print("   -> Processed ", i, " / ", #zerosL));

    gam = zerosL[i];
    rho = 0.5 + I*gam;
    
    zeta_rho = lfun(Z, rho);          \\ zeta(rho)
    Lprime_rho = lfun(L, rho, 1);     \\ L'(rho, chi_-4)
    
    listput(linha1, Strprintf("%.30g %.30g %.30g", gam, real(zeta_rho), imag(zeta_rho)));
    listput(linha2, Strprintf("%.30g %.30g %.30g", gam, real(Lprime_rho), imag(Lprime_rho)));
});

\\ Write the files (using Vec to convert List to vector)
write("zeta_zeros_L4_prec.txt", Str("# gam   Re(zeta(rho))   Im(zeta(rho))\n", strjoin(Vec(linha1), "\n")));
write("Lprime_zeros_L4_prec.txt", Str("# gam   Re(Lprime(rho))   Im(Lprime(rho))\n", strjoin(Vec(linha2), "\n")));
print("L(s,chi_-4) data written.");

\\ =========================================================================
\\ BLOCK 2: Zeros of Riemann zeta
\\ For these zeros:  zeta_Q(i)'(rho) = zeta'(rho) * L(rho)
\\ =========================================================================
linha3 = List();
linha4 = List();

print("Processing Riemann zeta zeros...");
for(i=1, #zerosZeta, {
    if(i % 500 == 0, print("   -> Processed ", i, " / ", #zerosZeta));

    gam = zerosZeta[i];
    rho = 0.5 + I*gam;
    
    zetaprime_rho = lfun(Z, rho, 1);  \\ zeta'(rho)
    L_rho = lfun(L, rho);             \\ L(rho, chi_-4)
    
    listput(linha3, Strprintf("%.30g %.30g %.30g", gam, real(zetaprime_rho), imag(zetaprime_rho)));
    listput(linha4, Strprintf("%.30g %.30g %.30g", gam, real(L_rho), imag(L_rho)));
});

write("zetaprime_zeros_zeta_prec.txt", Str("# gam   Re(zetaprime(rho))   Im(zetaprime(rho))\n", strjoin(Vec(linha3), "\n")));
write("L_zeros_zeta_prec.txt", Str("# gam   Re(L(rho))   Im(L(rho))\n", strjoin(Vec(linha4), "\n")));
print("Riemann zeta data written.");

print("=========================================================================");
print("All files generated successfully up to height ", Tmax);
print("=========================================================================");