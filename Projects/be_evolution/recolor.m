#!/usr/bin/octave -q
imread("~/flower.png");

mx=500;
tr=15;

for f = 0:120
  try
    fff=sprintf("%03d", f);
    in=["raw_image_" fff ".txt"];
    load(in);
    I=flipud(image');
    M=(I>=0);
    I=I.*M;
    W=max(I-tr,0)/(mx-tr);
    B=min(I,tr)/tr;
    imwrite(cat(3,1-M,W,B),["recolor_" fff ".png"]);
  catch
    printf("fail on %d\n",f);
  end
endfor;


