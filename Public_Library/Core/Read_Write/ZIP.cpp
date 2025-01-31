//#####################################################################
// Copyright 2009-2010, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_ZLIB_SUPPORT
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/READ_WRITE.h>
#include <Core/Read_Write/ZIP.h>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
namespace PhysBAM{
//#####################################################################
// class ZIP_FILE_HEADER
//#####################################################################
struct ZIP_FILE_HEADER
{
    unsigned short version;
    unsigned short flags;
    unsigned short compression_type;
    unsigned short stamp_date,stamp_time;
    unsigned int crc;
    unsigned int compressed_size,uncompressed_size;
    std::string filename;
    unsigned int header_offset; // local header offset

    ZIP_FILE_HEADER()
    {}

    ZIP_FILE_HEADER(const std::string& filename_input)
        :version(20),flags(0),compression_type(8),stamp_date(0),stamp_time(0),crc(0),
        compressed_size(0),uncompressed_size(0),filename(filename_input),header_offset(0)
    {}

    bool Read(std::istream& istream,const bool global)
    {unsigned int sig;
    unsigned short version,flags;
    // read and check for local/global magic
    if(global){
        Read_Primitive(istream,sig);
        if(sig!=0x02014b50){LOG::cerr<<"Did not find global header signature"<<std::endl;return false;}
        Read_Primitive(istream,version);}
    else{
        Read_Primitive(istream,sig);
        if(sig!=0x04034b50){LOG::cerr<<"Did not find local header signature"<<std::endl;return false;}}
    // Read rest of header
    Read_Primitive(istream,version);
    Read_Primitive(istream,flags);
    Read_Primitive(istream,compression_type);
    Read_Primitive(istream,stamp_date);
    Read_Primitive(istream,stamp_time);
    Read_Primitive(istream,crc);
    Read_Primitive(istream,compressed_size);
    Read_Primitive(istream,uncompressed_size);
    unsigned short filename_length,extra_length;
    Read_Primitive(istream,filename_length);
    Read_Primitive(istream,extra_length);
    unsigned short comment_length=0;
    if(global){
        Read_Primitive(istream,comment_length); // filecomment
        unsigned short disk_number_start,int_file_attrib;
        unsigned int ext_file_attrib;
        Read_Primitive(istream,disk_number_start); // disk# start
        Read_Primitive(istream,int_file_attrib); // internal file
        Read_Primitive(istream,ext_file_attrib); // ext final
        Read_Primitive(istream,header_offset);} // rel offset
    char* buf=new char[std::max(comment_length,std::max(filename_length,extra_length))];
    istream.read(buf,filename_length);
    buf[filename_length]=0;
    filename=std::string(buf);
    istream.read(buf,extra_length);
    if(global) istream.read(buf,comment_length);
    delete [] buf;
    return true;}

    void Write(std::ostream& ostream,const bool global) const
    {if(global){
        Write_Primitive(ostream,(unsigned int)0x02014b50); // header sig
        Write_Primitive(ostream,(unsigned short)00);} // version made by
    else Write_Primitive(ostream,(unsigned int)0x04034b50);
    Write_Primitive(ostream,version);
    Write_Primitive(ostream,flags);
    Write_Primitive(ostream,compression_type);
    Write_Primitive(ostream,stamp_date);
    Write_Primitive(ostream,stamp_time);
    Write_Primitive(ostream,crc);
    Write_Primitive(ostream,compressed_size);
    Write_Primitive(ostream,uncompressed_size);
    Write_Primitive(ostream,(unsigned short)filename.length());
    Write_Primitive(ostream,(unsigned short)0); // extra lengthx
    if(global){
        Write_Primitive(ostream,(unsigned short)0); // filecomment
        Write_Primitive(ostream,(unsigned short)0); // disk# start
        Write_Primitive(ostream,(unsigned short)0); // internal file
        Write_Primitive(ostream,(unsigned int)0); // ext final
        Write_Primitive(ostream,(unsigned int)header_offset);} // rel offset
    for(unsigned int i=0;i<filename.length();i++) Write_Primitive(ostream,filename.c_str()[i]);}
//#####################################################################
};
//#####################################################################
// class GZIP_FILE_HEADER
//#####################################################################
struct GZIP_FILE_HEADER
{
    unsigned char magic0,magic1; // magic should be 0x8b,0x1f
    unsigned char cm; // compression method 0x8 is gzip
    unsigned char flags; // flags
    unsigned int modtime; // 4 byte modification time
    unsigned char flags2; // secondary flags
    unsigned char os; // operating system 0xff for unknown
    unsigned short crc16; // crc check
    unsigned int crc32;

    GZIP_FILE_HEADER()
        :magic0(0),magic1(0),flags(0),modtime(0),flags2(0),os(0),crc16(0),crc32(0)
    {}

    bool Read(std::istream& istream)
    {Read_Primitive(istream,magic0);
    Read_Primitive(istream,magic1);
    if(magic0 != 0x1f || magic1 != 0x8b){LOG::cerr<<"gzip: did not find gzip magic 0x1f 0x8b"<<std::endl;return false;}
    Read_Primitive(istream,cm);
    if(cm!=8){LOG::cerr<<"gzip: compression method not 0x8"<<std::endl;return false;}
    Read_Primitive(istream,flags);
    Read_Primitive(istream,modtime);
    Read_Primitive(istream,flags2);
    Read_Primitive(istream,os);
    unsigned char dummyByte;
    // read flags if necessary
    if(flags&2){
        unsigned short flgExtraLen;
        Read_Primitive(istream,flgExtraLen);
        for(int k=0;k<flgExtraLen;k++) Read_Primitive(istream,dummyByte);}
    // read filename/comment if present
    int stringsToRead=((flags&8)?1:0) + ((flags&4)?1:0);
    for(int i=0;i<stringsToRead;i++) 
        do{Read_Primitive(istream,dummyByte);}while(dummyByte!=0 && istream);
    if(flags&1) Read_Primitive(istream,crc16);
    if(!istream) {LOG::cerr<<"gzip: got to end of file after only reading gzip header"<<std::endl;return false;}
    return true;}

    void Write(std::ostream& ostream)
    {magic0=0x1f;magic1=0x8b;cm=8;flags=0;os=0xff;
    Write_Primitive(ostream,magic0);Write_Primitive(ostream,magic1);
    Write_Primitive(ostream,cm);
    Write_Primitive(ostream,flags);
    Write_Primitive(ostream,modtime);
    Write_Primitive(ostream,flags2);
    Write_Primitive(ostream,os);}

//#####################################################################
};
//#####################################################################
// class ZIP_STREAMBUF_DECOMPRESS
//#####################################################################
class ZIP_STREAMBUF_DECOMPRESS:public std::streambuf
{
    static const unsigned int buffer_size=512;
    std::istream& istream;

    z_stream strm;
    unsigned char in[buffer_size],out[buffer_size];
    ZIP_FILE_HEADER header;
    GZIP_FILE_HEADER gzip_header;
    int total_read,total_uncompressed;
    bool part_of_zip_file;
    bool valid;
    bool compressed_data;

    static const unsigned short DEFLATE=8;
    static const unsigned short UNCOMPRESSED=0;
public:
    ZIP_STREAMBUF_DECOMPRESS(std::istream& stream,bool part_of_zip_file_input)
        :istream(stream),total_read(0),total_uncompressed(0),part_of_zip_file(part_of_zip_file_input),valid(true)
    {
        strm.zalloc=Z_NULL;strm.zfree=Z_NULL;strm.opaque=Z_NULL;strm.avail_in=0;strm.next_in=Z_NULL;
        setg((char*)in,(char*)in,(char*)in);
        setp(0,0);
        // skip the header
        if(part_of_zip_file){
            valid=header.Read(istream,false);
            if(header.compression_type==DEFLATE) compressed_data=true;
            else if(header.compression_type==UNCOMPRESSED) compressed_data=false;
            else{
                compressed_data=false;LOG::cerr<<"ZIP: got unrecognized compressed data (Supported deflate/uncompressed)"<<std::endl;
                valid=false;}}
        else{valid=gzip_header.Read(istream);compressed_data=true;}
        // initialize the inflate
        if(compressed_data && valid){
            int result=inflateInit2(&strm,-MAX_WBITS);
            if(result!=Z_OK){LOG::cerr<<"gzip: inflateInit2 did not return Z_OK"<<std::endl;valid=false;}}
    }

    virtual ~ZIP_STREAMBUF_DECOMPRESS()
    {if(compressed_data && valid) inflateEnd(&strm);
    if(!part_of_zip_file) delete &istream;}

    int process()
    {if(!valid) return -1;
    if(compressed_data){
        strm.avail_out=buffer_size-4;
        strm.next_out=(Bytef*)(out+4);
        while(strm.avail_out!=0){
            if(strm.avail_in==0){ // buffer empty, read some more from file
                istream.read((char*)in,part_of_zip_file?std::min((unsigned int)buffer_size,header.compressed_size-total_read):(unsigned int)buffer_size);
                strm.avail_in=istream.gcount();
                total_read+=strm.avail_in;
                strm.next_in=(Bytef*)in;}
            int ret=inflate(&strm,Z_NO_FLUSH); // decompress
            switch(ret){
                case Z_STREAM_ERROR: 
                    LOG::cerr<<"libz error Z_STREAM_ERROR"<<std::endl;
                    valid=false;return -1;
                case Z_NEED_DICT:
                case Z_DATA_ERROR:
                case Z_MEM_ERROR: 
                    LOG::cerr<<"gzip error "<<strm.msg<<std::endl;
                    valid=false;return -1;}
            if(ret==Z_STREAM_END) break;}
        int unzip_count=buffer_size-strm.avail_out-4;
        total_uncompressed+=unzip_count;
        return unzip_count;}
    else{ // uncompressed, so just read
        istream.read((char*)(out+4),std::min(buffer_size-4,header.uncompressed_size-total_read));
        int count=istream.gcount();
        total_read+=count;
        return count;}
    return 1;}

    virtual int underflow()
    {if(gptr() && (gptr()<egptr())) return traits_type::to_int_type(*gptr()); // if we already have data just use it
    int put_back_count=gptr()-eback();
    if(put_back_count>4) put_back_count=4;
    std::memmove(out+(4-put_back_count),gptr()-put_back_count,put_back_count);
    int num=process();
    setg((char*)(out+4-put_back_count),(char*)(out+4),(char*)(out+4+num));
    if(num<=0) return EOF;
    return traits_type::to_int_type(*gptr());}

    virtual int overflow(int c=EOF)
    {assert(false);return EOF;}

//#####################################################################
};

//#####################################################################
// class ZIP_STREAMBUF_COMPRESS
//#####################################################################
class ZIP_STREAMBUF_COMPRESS:public std::streambuf
{
    static const int buffer_size=512;
    std::ostream& ostream; // owned when header==0 (when not part of zip file)

    z_stream strm;
    unsigned char in[buffer_size],out[buffer_size];

    ZIP_FILE_HEADER* header;
    GZIP_FILE_HEADER gzip_header;
    unsigned int header_offset;
    unsigned int uncompressed_size;
    unsigned int crc;

    bool valid;

public:
    ZIP_STREAMBUF_COMPRESS(ZIP_FILE_HEADER* header,std::ostream& stream)
        :ostream(stream),header(header),valid(true)
    {
        strm.zalloc=Z_NULL;strm.zfree=Z_NULL;strm.opaque=Z_NULL;
        int ret=deflateInit2(&strm,Z_DEFAULT_COMPRESSION,Z_DEFLATED,-MAX_WBITS,8,Z_DEFAULT_STRATEGY); 
        if(ret != Z_OK){LOG::cerr<<"libz: failed to deflateInit"<<std::endl;valid=false;return;}
        setg(0,0,0);
        setp((char*)in,(char*)(in+buffer_size-4)); // we want to be 4 aligned
        // Write appropriate header
        if(header){header->header_offset=stream.tellp();header->Write(ostream,false);}
        else{header_offset=stream.tellp();gzip_header.Write(ostream);}
        uncompressed_size=crc=0;
    }

    virtual ~ZIP_STREAMBUF_COMPRESS()
    {if(valid){
        process(true);
        deflateEnd(&strm);
        if(header){
            std::ios::streampos final_position=ostream.tellp();
            header->uncompressed_size=uncompressed_size;
            header->crc=crc;
            ostream.seekp(header->header_offset);
            header->Write(ostream,false);
            ostream.seekp(final_position);}
        else{Write_Primitive(ostream,crc);Write_Primitive(ostream,uncompressed_size);}}
    if(!header) delete &ostream;}

protected:
    int process(bool flush)
    {if(!valid) return -1;
    strm.next_in=(Bytef*)pbase();
    strm.avail_in=pptr()-pbase();
    while(strm.avail_in!=0 || flush){
        strm.avail_out=buffer_size;
        strm.next_out=(Bytef*)out;
        int ret=deflate(&strm,flush?Z_FINISH:Z_NO_FLUSH);
        if(!(ret!=Z_BUF_ERROR && ret!=Z_STREAM_ERROR)){
            valid=false;
            LOG::cerr<<"gzip: gzip error "<<strm.msg<<std::endl;;
            return -1;}
        int generated_output=strm.next_out-(Bytef*)out;
        ostream.write((char*)out,generated_output);
        if(header) header->compressed_size+=generated_output;
        if(ret==Z_STREAM_END) break;}
    // update counts, crc's and buffers
    int consumed_input=pptr()-pbase();
    uncompressed_size+=consumed_input;
    crc=crc32(crc,(Bytef*)in,consumed_input);
    setp(pbase(),pbase()+buffer_size-4);return 1;}

    virtual int sync()
    {if(pptr() && pptr()>pbase()) return process(false);return 0;}

    virtual int underflow()
    {PHYSBAM_FATAL_ERROR("Attempt to read write only ostream");return 0;}

    virtual int overflow(int c=EOF)
    {if(c!=EOF){*pptr()=c;pbump(1);}
    if(process(false)==EOF) return EOF;
    return c;}

//#####################################################################
};
//#####################################################################
// Class ZIP_FILE_ISTREAM
//#####################################################################
// Class needed because istream cannot own its streambuf
class ZIP_FILE_ISTREAM:public std::istream
{
    ZIP_STREAMBUF_DECOMPRESS buf;
public:
    ZIP_FILE_ISTREAM(std::istream& istream,bool part_of_zip_file)
        :std::istream(&buf),buf(istream,part_of_zip_file)
    {}

    virtual ~ZIP_FILE_ISTREAM()
    {}

//#####################################################################
};
//#####################################################################
// Class ZIP_FILE_OSTREAM
//#####################################################################
// Class needed because ostream cannot own its streambuf
class ZIP_FILE_OSTREAM:public std::ostream
{
    ZIP_STREAMBUF_COMPRESS buf;
public:
    ZIP_FILE_OSTREAM(ZIP_FILE_HEADER* header,std::ostream& ostream)
        :std::ostream(&buf),buf(header,ostream)
    {}

    virtual ~ZIP_FILE_OSTREAM()
    {}

//#####################################################################
};
//#####################################################################
// Function ZIP_FILE_WRITER
//#####################################################################
ZIP_FILE_WRITER::
ZIP_FILE_WRITER(const std::string& filename)
{
    ostream.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(!ostream) throw FILESYSTEM_ERROR("ZIP: Invalid file handle");
}
//#####################################################################
// Function ZIP_FILE_WRITER
//#####################################################################
ZIP_FILE_WRITER::
~ZIP_FILE_WRITER()
{
    // Write all file headers
    std::ios::streampos final_position=ostream.tellp();
    for(int i=0;i<files.Size();i++){files(i)->Write(ostream,true);delete files(i);}
    std::ios::streampos central_end=ostream.tellp();
    // Write end of central
    Write_Primitive(ostream,(unsigned int)0x06054b50); // end of central
    Write_Primitive(ostream,(unsigned short)0); // this disk number
    Write_Primitive(ostream,(unsigned short)0); // this disk number
    Write_Primitive(ostream,(unsigned short)files.Size()); // one entry in center in this disk
    Write_Primitive(ostream,(unsigned short)files.Size()); // one entry in center 
    Write_Primitive(ostream,(unsigned int)(central_end-final_position)); // size of header
    Write_Primitive(ostream,(unsigned int)final_position); // offset to header
    Write_Primitive(ostream,(unsigned short)0); // zip comment
}
//#####################################################################
// Function ZIP_FILE_WRITER
//#####################################################################
std::ostream* ZIP_FILE_WRITER::
Add_File(const std::string& filename,const bool binary)
{
    files.Append(new ZIP_FILE_HEADER(filename));
    return new ZIP_FILE_OSTREAM(files.Last(),ostream);
}
//#####################################################################
// Function ZIP_FILE_READER
//#####################################################################
ZIP_FILE_READER::
ZIP_FILE_READER(const std::string& filename)
{
    istream.open(filename.c_str(),std::ios::in|std::ios::binary);
    if(!istream) throw FILESYSTEM_ERROR("ZIP: Invalid file handle");
    Find_And_Read_Central_Header();
}
//#####################################################################
// Function ZIP_FILE_READER
//#####################################################################
ZIP_FILE_READER::
~ZIP_FILE_READER()
{
    for(auto&v:filename_to_header){delete v.data;v.data=0;}
}
//#####################################################################
// Function Find_And_Read_Central_Header
//#####################################################################
bool ZIP_FILE_READER::
Find_And_Read_Central_Header()
{
    // Find the header
    // NOTE: this assumes the zip file header is the last thing written to file...
    istream.seekg(0,std::ios_base::end);
    std::ios::streampos end_position=istream.tellg();
    unsigned int max_comment_size=0xffff; // max size of header
    unsigned int read_size_before_comment=22;
    std::ios::streamoff read_start=max_comment_size+read_size_before_comment;
    if(read_start>end_position) read_start=end_position;
    istream.seekg(end_position-read_start);
    if(read_start<=0){LOG::cerr<<"ZIP: Invalid read buffer size"<<std::endl;return false;}
    char *buf=new char[read_start];
    istream.read(buf,read_start);
    int found=-1;
    for(unsigned int i=0;(int)i<read_start-3;i++){
        if(buf[i]==0x50 && buf[i+1]==0x4b && buf[i+2]==0x05 && buf[i+3]==0x06){found=i;break;}}
    delete [] buf;
    if(found==-1){LOG::cerr<<"ZIP: Failed to find zip header"<<std::endl;return false;}
    // seek to end of central header and read
    istream.seekg(end_position-(read_start-found));
    unsigned int word;
    unsigned short disk_number1,disk_number2,num_files,num_files_this_disk;
    Read_Primitive(istream,word); // end of central
    Read_Primitive(istream,disk_number1); // this disk number
    Read_Primitive(istream,disk_number2); // this disk number
    if(disk_number1!=disk_number2 || disk_number1!=0){
        LOG::cerr<<"ZIP: multiple disk zip files are not supported"<<std::endl;return false;}
    Read_Primitive(istream,num_files); // one entry in center in this disk
    Read_Primitive(istream,num_files_this_disk); // one entry in center 
    if(num_files != num_files_this_disk){
        LOG::cerr<<"ZIP: multi disk zip files are not supported"<<std::endl;return false;}
    unsigned int size_of_header,header_offset;
    Read_Primitive(istream,size_of_header); // size of header
    Read_Primitive(istream,header_offset); // offset to header
    // go to header and read all file headers
    istream.seekg(header_offset);
    for(int i=0;i<num_files;i++){
        ZIP_FILE_HEADER* header=new ZIP_FILE_HEADER;
        bool valid=header->Read(istream,true);
        if(valid) filename_to_header.Insert(header->filename,header);
        else delete header;}
    return true;
}
//#####################################################################
// Function Get_File
//#####################################################################
std::istream* ZIP_FILE_READER::Get_File(const std::string& filename,const bool binary)
{
    ZIP_FILE_HEADER** header=filename_to_header.Get_Pointer(filename);
    if(header){istream.seekg((*header)->header_offset);return new ZIP_FILE_ISTREAM(istream,true);}
    return 0;
}
//#####################################################################
// Function Get_File_List
//#####################################################################
void ZIP_FILE_READER::Get_File_List(ARRAY<std::string>& filenames) const
{
    filenames.Remove_All();
    filename_to_header.Get_Keys(filenames);
}
//#####################################################################
// Function Gzip_In
//#####################################################################
std::istream* 
Gzip_In(const std::string& filename,std::ios::openmode mode)
{
    std::ifstream* infile=new std::ifstream(filename.c_str(),mode | std::ios::binary);
    return new ZIP_FILE_ISTREAM(*infile,false);
}
//#####################################################################
// Function Gzip_Out
//#####################################################################
std::ostream* 
Gzip_Out(const std::string& filename,std::ios::openmode mode)
{
    std::ofstream* outfile=new std::ofstream(filename.c_str(),mode);
    return new ZIP_FILE_OSTREAM(0,*outfile);
}
//#####################################################################
}
#endif
