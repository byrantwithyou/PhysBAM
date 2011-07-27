#!/usr/bin/octave -q
load(argv(){1});
imwrite([argv(){1} ".png"], EI(:,:,1), EI(:,:,2), EI(:,:,3));
