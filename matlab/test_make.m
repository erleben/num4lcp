% Copyright 2011, Kenny Erleben (DIKU)
clear all
close all
clc

A = make_pd_matrix(50);
[x, b] = make_lcp(A,0.5);
tst = x'*(A*x+b)
figure(1)
spy(A);
title('Random PD matrix');
figure(2)
plot(sort(eig(A)));
title('Eigenvalues of random PD matrix');

A = make_psd_matrix(50,0.5);
[x, b] = make_lcp(A,0.5);
tst = x'*(A*x+b)
figure(3)
spy(A);
title('Random PSD matrix');
figure(4)
plot(sort(eig(A)));
title('Eigenvalues of random PSD matrix');

A = make_blocked_matrix(3,30,0.5);
[x, b] = make_lcp(A,0.5);
tst = x'*(A*x+b)
figure(5)
spy(A);
title('Random blocked matrix');
figure(6)
plot(sort(eig(A)));
title('Eigenvalues of random blocked matrix');

A = make_banded_matrix(3,100);
[x, b] = make_lcp(A,0.5);
tst = x'*(A*x+b)
figure(7)
spy(A);
title('Random banded matrix');
figure(8)
plot(sort(eig(A)));
title('Eigenvalues of random banded matrix');

A = make_contact_matrix(10);
[x, b] = make_lcp(A,0.5);
tst = x'*(A*x+b)
figure(9)
spy(A);
title('Random contact matrix');
print('-f9','-depsc2', 'output/contact_matrix');
figure(10)
plot(sort(eig(A)));
title('Eigenvalues of random contact matrix');

A = make_fluid_matrix(4);
[x, b] = make_lcp(A,0.5);
tst = x'*(A*x+b)
figure(11)
spy(A);
title('Random fluid matrix');
print('-f11','-depsc2', 'output/fluid_matrix');
figure(12)
plot(sort(eig(A)));
title('Eigenvalues of random fluid matrix');
