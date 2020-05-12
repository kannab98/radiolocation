clc, clear
Fd = 1e6;
t = 0:1/Fd:10e-6;
fn = 100e3;
W = exp(1i*2*pi*fn*t).';
figure(1), plot(t, real(W))
f = 0:10e2:200e3;
n = length(f);
out=zeros(1,n);
for ii=1:n
    signal = exp(1i*2*pi*f(1,ii)*t).';
    out(1,ii) = W'*signal;
end
figure(2), plot(f, 20*log10(abs(out)))