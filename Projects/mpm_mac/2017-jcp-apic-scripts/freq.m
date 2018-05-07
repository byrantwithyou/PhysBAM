function freq = freq(name,start)
a = importdata(strcat(name,'-time'));
b = importdata(strcat(name,'-vort'));
a = a(start:end);
b = b(start:end);
[v,t] = resample(b,a);
f = 1/(t(2)-t(1));

psdest = psd(spectrum.periodogram,v,'Fs',f,'NFFT',length(v));
[~,I] = max(psdest.Data);
fprintf('Maximum occurs at %d Hz.\n',psdest.Frequencies(I));
fprintf('Period: %f.\n',1/psdest.Frequencies(I));

