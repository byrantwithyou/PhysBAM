#include <stdio.h>
#include <stdlib.h>


void Check_Allocatable_Memory()
{
    const int CHUNK_SIZE=1024*1024; // default 1MB
    const int MAX_NUMBER_OF_CHUNKS=32000; // default 32 GB
    char** temp_buffer=new char*[MAX_NUMBER_OF_CHUNKS];
    for(int i=0;i<MAX_NUMBER_OF_CHUNKS;i++){
        temp_buffer[i]=new char[CHUNK_SIZE];
        printf("Allocated %d Megabytes.\n",i);
    }
}

int main(int argc,char* argv[])
{
    printf(
           "Performance Tester. Copyright (c) 2005, Frank Losasso\n"\
           "\n"\
           "  1) Check Allocatable Memory Amount\n");
    int choice=0;
    scanf("%d",&choice);
    
    switch(choice){
      case 1:
        Check_Allocatable_Memory();break;
      default:
        printf("Invalid Choice\n");exit(1);
    }
    exit(0);
}
