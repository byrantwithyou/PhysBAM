#!/usr/bin/python
import Image
import ImageDraw
import ImageFont
import shutil
import rangeutils
import os
import shutil

def slide(filename,size,lines):
    image_size=640,480
    img=Image.new('RGB',image_size)
    draw=ImageDraw.Draw(img)

    fontname='/u/aselle/shared/etc/VeraBd.ttf'
    font=ImageFont.truetype(fontname,size)

    # compute line length
    offsets=[]
    widths=[]
    fonts=[]
    colors=[]
    height=0
    default_color='#ffeedd'
    for line in lines:
        color=default_color
        if type(line)==tuple:
            fonts.append(ImageFont.truetype(fontname,line[0]))
            if len(line)>2:
                color=line[2]
            line=line[1]
        else:
            fonts.append(font)
            
        w,h=draw.textsize(line,font=fonts[-1])
        offsets.append(height)
        widths.append(w)
        colors.append(color)
        height+=int(h+.10*h)

    height_start=image_size[1]/2-height/2
    for index in range(len(offsets)):
        if type(lines[index])==tuple:
            line=lines[index][1]
        else:
            line=lines[index]
        draw.text((image_size[0]/2-widths[index]/2,height_start+offsets[index]),line,font=fonts[index],fill=colors[index])

    img.save(filename)

class VIDEO:
    def __init__(self,directory):
        self.directory=directory
        os.mkdir(self.directory)
        self.frame=0
        self.fps=30

    def frame_file(self,count):
        return os.path.join(self.directory,"video.%06d.png"%self.frame)

    def add_frame(self,image_filename,count=1):
        #print "adding frame %s with %d copies (index=%d)"%(image_filename,count,self.frame),
        if count>1: print "Adding %-90s %3.1f s"%(image_filename,count*self.fps)
        for i in range(count):
            #print "Copying %s to %s"%(image_filename,self.frame_file(self.frame))
            #os.symlink(image_filename,self.frame_file(self.frame))
            shutil.copy(image_filename,self.frame_file(self.frame))
            self.frame+=1
        #print self.frame

    def add_directory(self,directory,start=-1,end=-1):
        files=os.listdir(directory)
        files=filter(lambda x:x.endswith('.png'),files)
        files.sort()
        files=map(lambda x:os.path.join(directory,x),files)
        for file in files:
                    self.add_frame(file)
        
    def make_movie(self,filename):
#         cmd="ffmpeg -y -r 30 -i %s/video.%%06d.png -vcodec mjpeg -b 15000k -r %d %s_mjpeg.avi"%(self.directory,self.fps,filename)
#         os.system(cmd)
#         cmd="ffmpeg -y -r 30 -i %s/video.%%06d.png -vcodec xvid -b 3430k -r %d %s_divx.avi"%(self.directory,self.fps,filename)
#         os.system(cmd)
#         cmd="ffmpeg -y -r 30 -i %s/video.%%06d.png -vcodec mpeg4 -b 3430k -r %d %s.mov"%(self.directory,self.fps,filename)
#         os.system(cmd)
        cmd="ffmpeg -y -r 30 -i %s/video.%%06d.png -vcodec mpeg4 -b 1215k -r %d %s_low.mov"%(self.directory,self.fps,filename)
        os.system(cmd)
#        cmd="ffmpeg -y -r 30 -i %s/video.%%06d.png -vcodec rawvideo -b 3430k -r %d %s_raw.mov"%(self.directory,self.fps,filename)
#        os.system(cmd)

