//#####################################################################
// Copyright 2006-2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIDEO
//#####################################################################
#define __STDC_CONSTANT_MACROS
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <cstring>
#include "VIDEO.h"
#include <boost/format.hpp>
using namespace PhysBAM;

#ifdef USE_FFMPEG
extern "C"{
#include <ffmpeg/avcodec.h>
#if LIBAVCODEC_VERSION_INT>3344384
#include <stdint.h>
#define FFMPEG_USE_SWSCALE //TODO: Remove when windows version updated
#include <ffmpeg/swscale.h>
#endif
#include <ffmpeg/avformat.h>
}

template<class T>
class FFMPEG_VIDEO_READER:public VIDEO_READER<T>
{
private:
    AVFormatContext *format_context;
    AVCodecContext *codec_context;
    AVStream* stream;
    AVCodec* codec;
    AVFrame* yuv_frame;
    AVFrame* rgb_frame;
    unsigned char* set_image;
    int stream_index;
    int current_frame;
    bool failed;
#ifdef FFMPEG_USE_SWSCALE
    struct SwsContext *img_convert_context;
#endif

public:
    FFMPEG_VIDEO_READER(const std::string& filename)
        :format_context(0),codec_context(0),stream(0),codec(0),yuv_frame(0),rgb_frame(0),set_image(0),stream_index(-1),current_frame(0),failed(false)
#ifdef FFMPEG_USE_SWSCALE
        ,img_convert_context(0)
#endif
    {
        av_register_all();
        if(av_open_input_file(&format_context,filename.c_str(),NULL,0,NULL)!=0){failed=true;return;}
        if(av_find_stream_info(format_context)<0){failed=true;return;}
        dump_format(format_context,0,filename.c_str(),false);
        for(unsigned int i=0;i<(unsigned int)format_context->nb_streams;i++) if(format_context->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO){stream_index=i;break;}
        if(stream_index==-1){failed=true;return;}
        stream=format_context->streams[stream_index];
        codec_context=format_context->streams[stream_index]->codec;
        codec=avcodec_find_decoder(codec_context->codec_id);
        if(!codec){failed=true;return;}
        //if(codec->capabilities & CODEC_CAP_TRUNCATED) codec_context->flags|=CODEC_FLAG_TRUNCATED;
        if(avcodec_open(codec_context,codec)<0){failed=true;return;}
        yuv_frame=avcodec_alloc_frame();
        rgb_frame=avcodec_alloc_frame();

        if(codec_context->time_base.den>1000 && codec_context->time_base.num==1) codec_context->time_base.num=1000;
        if(!stream->nb_frames) stream->nb_frames=(int)stream->duration;
    }

    ~FFMPEG_VIDEO_READER()
    {
        if(yuv_frame) av_free(yuv_frame);
        if(rgb_frame) av_free(rgb_frame);
        if(codec_context) avcodec_close(codec_context);
        if(format_context) av_close_input_file(format_context);
#ifdef FFMPEG_USE_SWSCALE
        if(img_convert_context) sws_freeContext(img_convert_context);
#endif
    }

    bool Valid() const
    {return !failed;}

    int Width() const
    {return codec_context->width;}

    int Height() const
    {return codec_context->height;}

    int Frame_Rate() const
    {return (int)((float)stream->time_base.den/codec_context->time_base.num);}

    int Frame_Count() const
    {return (int)stream->nb_frames;}

    void Frame(const int frame,unsigned char* image)
    {AVPacket packet;
    if(frame == 1){av_seek_frame(format_context,stream_index,0,AVSEEK_FLAG_ANY);}
    else if(frame != current_frame+1){
        av_seek_frame(format_context,stream_index,frame-1,AVSEEK_FLAG_ANY); }
    if(image != set_image){set_image=image;avpicture_fill((AVPicture*)rgb_frame,image,PIX_FMT_RGB24,Width(),Height());}
    while(av_read_frame(format_context,&packet)>=0){
        if(packet.stream_index==stream_index){
            int frame_finished=0;
            avcodec_decode_video(codec_context,yuv_frame,&frame_finished,packet.data,packet.size);
            av_free_packet(&packet);
            if(frame_finished){
                current_frame=frame;
#ifdef FFMPEG_USE_SWSCALE
                if(!img_convert_context){
                    img_convert_context=sws_getContext(Width(),Height(),codec_context->pix_fmt,Width(),Height(),PIX_FMT_RGB24,SWS_POINT,NULL,NULL,NULL);
                    if(!img_convert_context) PHYSBAM_FATAL_ERROR("ffmpeg: Failed to allocate scaler context");}
                sws_scale(img_convert_context,yuv_frame->data,yuv_frame->linesize,0,Height(),rgb_frame->data,rgb_frame->linesize);
#else
                img_convert((AVPicture*)rgb_frame,PIX_FMT_RGB24,(AVPicture*)yuv_frame,codec_context->pix_fmt,codec_context->width,codec_context->height);
#endif
                return;}}
        else av_free_packet(&packet);}}

    void Frame(const int frame,ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image)
    {int bytes_per_frame=avpicture_get_size(PIX_FMT_RGB24,Width(),Height());
    unsigned char* frame_buffer=new unsigned char[bytes_per_frame];
    Frame(frame,frame_buffer);
    image.Resize(1,Width(),1,Height(),false,false);
    image.Fill(VECTOR<T,3>());
    unsigned char* img=frame_buffer;
    for(int j=1;j<=Height();j++) for(int i=1;i<=Width();i++){
        image(i,j)=VECTOR<T,3>(IMAGE<T>::Byte_Color_To_Scalar_Color(img[0]),IMAGE<T>::Byte_Color_To_Scalar_Color(img[1]),IMAGE<T>::Byte_Color_To_Scalar_Color(img[2]));img+=3;}

    delete[] frame_buffer;}

//#####################################################################
};
template class FFMPEG_VIDEO_READER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FFMPEG_VIDEO_READER<double>; 
#endif
#endif

//#####################################################################
// Class IMAGE_SEQUENCE
//#####################################################################
template<class T>
class IMAGE_SEQUENCE_READER:public VIDEO_READER<T>
{
private:
    std::string filename;
    int min_frame,max_frame;
    int width,height;
    ARRAY<VECTOR<T,3>,VECTOR<int,2> > temporary_image;
    ARRAY<unsigned char*> cache;

public:
    IMAGE_SEQUENCE_READER(const std::string& filename_input,const int min_frame_input,const int max_frame_input)
        :filename(filename_input),min_frame(min_frame_input),max_frame(max_frame_input)
    {
        Frame(1,temporary_image);
        width=temporary_image.counts.x;height=temporary_image.counts.y;
        cache.Resize(Frame_Count());
    }

    ~IMAGE_SEQUENCE_READER()
    {
        for(int i=1;i<=cache.m;i++) if(cache(i)) delete[] cache(i);
    }

    bool Valid() const
    {return true;}

    int Width() const
    {return width;}

    int Height() const
    {return height;}

    int Frame_Rate() const
    {return 24;}

    int Frame_Count() const
    {return max_frame-min_frame+1;}

    void Frame(const int frame,unsigned char* image)
    {int buffer_size=Width()*Height()*3;
    if(cache(frame)) memcpy(image,cache(frame),buffer_size);
    else{
        Frame(frame,temporary_image);
        unsigned char* buffer=image;
        for(int j=height;j>=1;j--) for(int i=1;i<=width;i++){
            buffer[0]=IMAGE<T>::Scalar_Color_To_Byte_Color(temporary_image(i,j).x);
            buffer[1]=IMAGE<T>::Scalar_Color_To_Byte_Color(temporary_image(i,j).y);
            buffer[2]=IMAGE<T>::Scalar_Color_To_Byte_Color(temporary_image(i,j).z);
            buffer+=3;}
        cache(frame)=new unsigned char[buffer_size];
        memcpy(cache(frame),image,buffer_size);}}
    

    void Frame(const int frame,ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image)
    {IMAGE<T>::Read(STRING_UTILITIES::string_sprintf(filename.c_str(),(frame+min_frame-1)),image);}

    //#####################################################################
};
template class IMAGE_SEQUENCE_READER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMAGE_SEQUENCE_READER<double>; 
#endif
//#####################################################################
// Function Read
//#####################################################################
template<class T> VIDEO_READER<T>* VIDEO_READER<T>::
Read(const std::string& filename)
{
    VIDEO_READER<T>* video=0;
#ifdef USE_FFMPEG
    video=new FFMPEG_VIDEO_READER<T>(filename);
#endif
    if(video && !video->Valid()){delete video;video=0;}
    return video;
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> VIDEO_READER<T>* VIDEO_READER<T>::
Read(const std::string& format_string,const int min_frame,const int max_frame)
{
    VIDEO_READER<T>* video=new IMAGE_SEQUENCE_READER<T>(format_string,min_frame,max_frame);
    if(!video->Valid()){delete video;video=0;}
    return video;
}
//#####################################################################


//#####################################################################
// Class FFMPEG_VIDEO_WRITER
//#####################################################################
#ifdef USE_FFMPEG
template<class T>
class FFMPEG_VIDEO_WRITER:public VIDEO_WRITER<T>
{
private:
    std::string filename;
    std::string codec_str;
    int bitrate;
    int frame_rate;
    int width,height;

    AVFormatContext* format;
    AVStream* stream;
    AVCodec* codec;
    AVFrame *picture_rgb,*picture_fmt;
    uint8_t *video_outbuf;
    int video_outbuf_size;
#ifdef FFMPEG_USE_SWSCALE
    struct SwsContext *img_convert_context;
#endif

public:
    static AVFrame* alloc_picture(int pix_fmt, int width, int height)
    {
        AVFrame *picture;
        uint8_t *picture_buf;
        int size;
        
        picture = avcodec_alloc_frame();
        if (!picture)
            return NULL;
        size = avpicture_get_size(pix_fmt, width, height);
        picture_buf = (uint8_t*) malloc(size);
        if (!picture_buf) {
            av_free(picture);
            return NULL;
        }
        avpicture_fill((AVPicture *)picture, picture_buf, 
                       pix_fmt, width, height);
        return picture;
    }

    FFMPEG_VIDEO_WRITER(const std::string& filename,const std::string& codec_str,const int bitrate,const int frame_rate,const int width,
        const int height)
        :filename(filename),codec_str(codec_str),bitrate(bitrate),frame_rate(frame_rate),width(width),height(height),format(0),stream(0),codec(0),picture_rgb(0),picture_fmt(0),
        video_outbuf(0),video_outbuf_size(0)
#ifdef FFMPEG_USE_SWSCALE
        ,img_convert_context(0)
#endif
    {
        av_register_all();

        //#############################################
        // Setup Format, Codec and parameters
        //#############################################
#if LIBAVFORMAT_VERSION_MAJOR<52
        format=av_alloc_format_context();
#else
        format=avformat_alloc_context();
#endif
        // output format and filename
        {AVOutputFormat* file_oformat=guess_format(NULL,filename.c_str(),NULL);
        if(!file_oformat) PHYSBAM_FATAL_ERROR("Can't find output format");
        format->oformat=file_oformat;}
        strncpy(format->filename,filename.c_str(),sizeof(format->filename));format->filename[sizeof(format->filename)-1]=0; // strncpy is broken
        format->timestamp=0;
        
        // video stream
        stream=av_new_stream(format,format->nb_streams);
        if(!stream) PHYSBAM_FATAL_ERROR("Could not allocate stream");
        avcodec_get_context_defaults(stream->codec);
        
        // lookup codec
        AVCodecContext* codec_context=stream->codec;
        codec_context->codec_id=av_guess_codec(format->oformat,NULL,format->filename,NULL,CODEC_TYPE_VIDEO);
        AVCodec* codec=avcodec_find_encoder_by_name(codec_str.c_str());
        if(codec==NULL) PHYSBAM_FATAL_ERROR(std::string("Unknown codec ")+codec_str);
        codec_context->codec_id=codec->id;
        codec_context->codec_type=CODEC_TYPE_VIDEO;
        // Frame rate & bit rate
        codec_context->time_base.den=frame_rate;
        codec_context->time_base.num=1;
        if(codec && codec->supported_framerates) PHYSBAM_NOT_IMPLEMENTED("Not implemented limited frame rates");
        codec_context->bit_rate=bitrate*1<<10;
        
        // size video
        codec_context->width=width;codec_context->height=height;
        codec_context->pix_fmt=PIX_FMT_YUV420P;
        if(codec && codec->pix_fmts){
            const enum PixelFormat* p=codec->pix_fmts;
            for(;*p!=-1;p++) if(*p==codec_context->pix_fmt) break;
            if(*p==-1){
                LOG::cerr<<"unsupported pixel format overriding"<<std::endl;
                codec_context->pix_fmt=codec->pix_fmts[0];}}

        if(av_set_parameters(format,NULL)<0) PHYSBAM_FATAL_ERROR("Invalid output format parameters");

        // open codec
        codec=avcodec_find_encoder(codec_context->codec_id);
        if(avcodec_open(codec_context,codec)<0) PHYSBAM_FATAL_ERROR("could not open codec");

        //#############################################
        // Make buffers
        //#############################################
        if(!(format->oformat->flags & AVFMT_RAWPICTURE)){
            video_outbuf_size=200000;
            video_outbuf=(uint8_t*) malloc(video_outbuf_size);}

        picture_rgb=alloc_picture(PIX_FMT_RGB24,width,height);
        picture_fmt=alloc_picture(codec_context->pix_fmt,width,height);
        if(picture_rgb==0 || picture_fmt==0) PHYSBAM_FATAL_ERROR("Failed to allocate picture buffers");
            

        //#############################################
        // Start writing
        //#############################################
        // open video file
        if (!(format->oformat->flags & AVFMT_NOFILE)){
            if (url_fopen(&format->pb, filename.c_str(), URL_WRONLY) < 0) PHYSBAM_FATAL_ERROR("could not open output file");}

        // write header
        dump_format(format,0,filename.c_str(),1);
        av_write_header(format);

    }

    virtual ~FFMPEG_VIDEO_WRITER()
    {
        // temporary buffers
        avcodec_close(stream->codec);
        av_free(picture_fmt->data[0]);
        av_free(picture_fmt);
        av_free(picture_rgb->data[0]);
        av_free(picture_rgb);
        av_free(video_outbuf);
        
#ifdef FFMPEG_USE_SWSCALE
        // free the scaler
        if(img_convert_context) sws_freeContext(img_convert_context);
#endif
        // write end
        av_write_trailer(format);

        // free streams
        for(unsigned int i=0;i<(unsigned int)format->nb_streams;i++) av_freep(&format->streams[i]);
#if LIBAVFORMAT_VERSION_MAJOR < 52
        if(!(format->oformat->flags & AVFMT_NOFILE)) url_fclose(&format->pb);
#else
        if(!(format->oformat->flags & AVFMT_NOFILE)) url_fclose(format->pb);
#endif
        
        // free format
        av_free(format);
    }

    int Frame(ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image)
    {
        int ret;

        unsigned char* buffer=picture_rgb->data[0];
        for(int j=height;j>=1;j--) for(int i=1;i<=width;i++){
            buffer[0]=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j).x);
            buffer[1]=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j).y);
            buffer[2]=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j).z);
            buffer+=3;}
#ifdef FFMPEG_USE_SWSCALE
        if(!img_convert_context){
            img_convert_context=sws_getContext(width,height,PIX_FMT_RGB24,width,height,stream->codec->pix_fmt,SWS_POINT,NULL,NULL,NULL);
            if(!img_convert_context) PHYSBAM_FATAL_ERROR("ffmpeg: Failed to allocate scaler context");}
        sws_scale(img_convert_context,picture_rgb->data,picture_rgb->linesize,0,height,picture_fmt->data,picture_fmt->linesize);
#else
        img_convert((AVPicture*)picture_fmt,stream->codec->pix_fmt,(AVPicture*)picture_rgb,PIX_FMT_RGB24,width,height);
#endif
        int out_size=avcodec_encode_video(stream->codec,video_outbuf,video_outbuf_size,picture_fmt);
        if(out_size!=0){
            AVPacket pkt;
            av_init_packet(&pkt);
            pkt.pts=stream->codec->coded_frame->pts;
            if(stream->codec->coded_frame->key_frame) pkt.flags |= PKT_FLAG_KEY;
            pkt.stream_index=stream->index;
            pkt.data=video_outbuf;
            pkt.size=out_size;
            ret=av_write_frame(format,&pkt);}
        else ret=0;
        return ret;
    }

    bool Valid() const
    {return true;}

};
template class FFMPEG_VIDEO_WRITER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FFMPEG_VIDEO_WRITER<double>; 
#endif
#endif

//#####################################################################
// Function Write
//#####################################################################
template<class T> VIDEO_WRITER<T>* VIDEO_WRITER<T>::
Write(const std::string& filename,const std::string& codec,const int bitrate,const int frame_rate,const int width,const int height)
{

#ifdef USE_FFMPEG
    VIDEO_WRITER<T>* video=new FFMPEG_VIDEO_WRITER<T>(filename,codec,bitrate,frame_rate,width,height);
    if(!video->Valid()){delete video;video=0;}
    return video;
#else    
    PHYSBAM_NOT_IMPLEMENTED();
#endif
}
//#####################################################################
// Function Enabled
//#####################################################################
template<class T> bool VIDEO_WRITER<T>::
Enabled()
{
#ifdef USE_FFMPEG
    return true;
#else
    return false;
#endif
}
//#####################################################################

template class VIDEO_READER<float>;
template class VIDEO_WRITER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VIDEO_READER<double>;
template class VIDEO_WRITER<double>;
#endif
