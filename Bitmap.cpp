//File: Bitmap.cpp
//Written by:     Mark Bernard
//on GameDev.net: Captain Jester
//e-mail: mark.bernard@rogers.com
//Please feel free to use and abuse this code as much
//as you like.  But, please give me some credit for
//starting you off on the right track.
//
//The file Bitmap.h goes along with this file
//
#include "Bitmap.h"
using namespace std;

//basic constructor
Bitmap::Bitmap(){
    reset();
}

//constructor loads the bitmap when it is created
Bitmap::Bitmap(char *file){
    reset();
    loadBMP(file);
}

//destructor
Bitmap::~Bitmap(){
    if(colours!=0) {
        delete[] colours;
    }
    if(data!=0) {
        delete[] data;
    }
}

//load a bitmap from a file and represent it correctly
//in memory
bool Bitmap::loadBMP(char *file) {
    FILE *in;                  //file stream for reading
    char *tempData;       //temp storage for image data
    int numColours;            //total available colours

    //bitmap is not loaded yet
    loaded=false;
    //make sure memory is not lost
    if(colours!=0) {
        delete[] colours;
    }
    if(data!=0) {
        delete[] data;
    }

    //open the file for reading in binary mode
    in=fopen(file,"rb");

    //if the file does not exist return in error
    if(in==NULL) {
        error="File not found";
        fclose(in);
        return false;
    }

    //read in the entire BITMAPFILEHEADER
    fread(&bmfh,sizeof(BitmapFileHeader),1,in);
	cout << "sizeof(BitmapFileHeader)=" << sizeof(BitmapFileHeader) << endl;
    //check for the magic number that says this is a bitmap
    if(bmfh.bfType!=BITMAP_MAGIC_NUMBER) {
        error="File is not in DIB format";
        fclose(in);
        return false;
    }

    //read in the entire BITMAPINFOHEADER
    fread(&bmih,sizeof(BitmapInfoHeader),1,in);
	cout << "sizeof(BitmapInfoHeader)=" << sizeof(BitmapInfoHeader) << endl;

    //save the width, height and bits per pixel for external use
    width=bmih.biWidth;
    height=bmih.biHeight;
    bpp=bmih.biBitCount;
	cout << "biBitCount      =" << bmih.biBitCount << endl;
	cout << "biClrImportant  =" << bmih.biClrImportant << endl;
	cout << "biClrUsed       =" << bmih.biClrUsed << endl;
	cout << "biCompression   =" << bmih.biCompression << endl;
	cout << "biHeight        =" << bmih.biHeight << endl;
	cout << "biPlanes        =" << bmih.biPlanes << endl;
	cout << "biSize          =" << bmih.biSize << endl;
	cout << "biSizeImage     =" << bmih.biSizeImage << endl;
	cout << "biWidth         =" << bmih.biWidth << endl;
	cout << "biXPelsPerMeter =" << bmih.biXPelsPerMeter << endl;
	cout << "biYPelsPerMeter =" << bmih.biYPelsPerMeter << endl;

    //calculate the size of the image data with padding
    dataSize=(width*height*(unsigned int)(bmih.biBitCount/8.0));

    //calculate the number of available colours
    numColours=1<<bmih.biBitCount;

    //if the bitmap is not 8 bits per pixel or more
    //return in error
    if(bpp<8) {
        error="File is not 8 or 24 bits per pixel";
        fclose(in);
        return false;
    }

    //load the palette for 8 bits per pixel
    if(bpp==8) {
    	colours=new RGBQuad[numColours];
    	fread(colours,sizeof(RGBQuad),numColours,in);
    }

    //set up the temporary buffer for the image data
    tempData=new char[dataSize];

    //exit if there is not enough memory
    if(tempData==NULL) {
        error="Not enough memory to allocate a temporary buffer";
        fclose(in);
        return false;
    }

    //read in the entire image
    fread(tempData,sizeof(char),dataSize,in);

    //close the file now that we have all the info
    fclose(in);

    //calculate the witdh of the final image in bytes
    byteWidth=padWidth=(int)((float)width*(float)bpp/8.0);

    //adjust the width for padding as necessary
    while(padWidth%4!=0) {
        padWidth++;
    }

    //change format from GBR to RGB
    if(bpp==8) {
    	loaded=convert8(tempData);
   	}
    else if(bpp==24) {
    	loaded=convert24(tempData);
   	}

    //clean up memory
    delete[] tempData;

    //bitmap is now loaded
    error="Bitmap loaded";

    //return success
    return loaded;
}

//function to set the inital values
void Bitmap::reset(void) {
    loaded=false;
    colours=0;
    data=0;
    error="";
}

bool Bitmap::convert24(char* tempData) {
	int offset,diff;

	diff=width*height*RGB_BYTE_SIZE;
    //allocate the buffer for the final image data
    data=new char[diff];

    //exit if there is not enough memory
    if(data==NULL) {
        error="Not enough memory to allocate an image buffer";
        delete[] data;
        return false;
    }

    if(height>0) {
        offset=padWidth-byteWidth;
        //count backwards so you start at the front of the image
        for(int i=0;i<dataSize;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+i+2)=*(tempData+i);
            *(data+i+1)=*(tempData+i+1);
            *(data+i)=*(tempData+i+2);
        }
    }

    //image parser for a forward image
    else {
        offset=padWidth-byteWidth;
        int j=dataSize-3;
        //count backwards so you start at the front of the image
		//here you can start from the back of the file or the front,
		//after the header  The only problem is that some programs
		//will pad not only the data, but also the file size to
		//be divisible by 4 bytes.
        for(int i=0;i<dataSize;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+j+2)=*(tempData+i);
            *(data+j+1)=*(tempData+i+1);
            *(data+j)=*(tempData+i+2);
            j-=3;
        }
    }

    return true;
}

bool Bitmap::convert8(char* tempData) {
	int offset,diff;

	diff=width*height*RGB_BYTE_SIZE;
    //allocate the buffer for the final image data
    data=new char[diff];

    //exit if there is not enough memory
    if(data==NULL) {
        error="Not enough memory to allocate an image buffer";
        delete[] data;
        return false;
    }

    if(height>0) {
        offset=padWidth-byteWidth;
        int j=0;
        //count backwards so you start at the front of the image
        for(int i=0;i<dataSize*RGB_BYTE_SIZE;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+i)=colours[*(tempData+j)].rgbRed;
            *(data+i+1)=colours[*(tempData+j)].rgbGreen;
            *(data+i+2)=colours[*(tempData+j)].rgbBlue;
            j++;
        }
    }

    //image parser for a forward image
    else {
        offset=padWidth-byteWidth;
        int j=dataSize-1;
        //count backwards so you start at the front of the image
        for(int i=0;i<dataSize*RGB_BYTE_SIZE;i+=3) {
            //jump over the padding at the start of a new line
            if((i+1)%padWidth==0) {
                i+=offset;
            }
            //transfer the data
            *(data+i)=colours[*(tempData+j)].rgbRed;
            *(data+i+1)=colours[*(tempData+j)].rgbGreen;
            *(data+i+2)=colours[*(tempData+j)].rgbBlue;
            j--;
        }
    }

    return true;
}

bool Bitmap::save8BMP(char * fileName,char * image, int width, int height)
{
	reset();
	saved = false;
	this->width = width;
	this->height = height;

	byteWidth = width;
	padWidth = width;
	while(padWidth%4!=0) padWidth++;


	//if(sizeof(&image))
    ofstream bmpfile(fileName, ios::out | ios::binary | ios::trunc);
    if(!bmpfile.is_open())
    {
        cout << "ERROR: FILE COULD NOT BE OPENED" << endl;
        return false;
    }
    CreateHeader(bmpfile);
    CreatePalette(bmpfile);
    CreatePixel(bmpfile, image);
    bmpfile.close();
    saved = true;
    return true;
}

void Bitmap::CreateHeader(ofstream &bfile)
{
    //int width = width;                        //will be coded with 4 bytes
    //int height = height;                       //... 4 bytes

	bmfh.bfType = 19778;                       //... 2 bytes
    bmfh.bfSize = (width*height*1)+40+14+1024;  //... 4 bytes //Replace 54 by 40
    bmfh.bfReserved1 = 0;                      //... 2 bytes
    bmfh.bfReserved2 = 0;                      //... 2 bytes
    bmfh.bfOffBits = 40+14+1024;                //... 4 bytes

    bmih.biSize = 40;                         //... 4 bytes //replace 54 by 40
    bmih.biWidth = width;						// 4 bytes
    bmih.biHeight = height;						// 4 bytes
    bmih.biPlanes = 1;                        //... 2 bytes
    bmih.biBitCount = 8;                           //... 2 bytes
    bmih.biCompression = 0;                    //... 4 bytes
    bmih.biSizeImage = (width*height*1);             //... 4 bytes
    bmih.biXPelsPerMeter = 2835;                          //... 4 bytes	0 default
    bmih.biYPelsPerMeter = 2835;                          //... 4 bytes	0 default
    bmih.biClrUsed = 256;                       //... 4 bytes
    bmih.biClrImportant = 256;                         //... 4 bytes 0 default

	// BMP FILE HEADER
    bfile.write((char*)&bmfh, sizeof(BitmapFileHeader));
    // BITMAP INFORMATION (DIB HEADER)
    bfile.write((char*)&bmih, sizeof(BitmapInfoHeader));
}

void Bitmap::CreatePalette(ofstream &bfile)
{
    int TabPalette[1024];
    int i=0;
    for (int x=0;x<256;x++)
    {
        TabPalette[i]=x;
        i++;
        TabPalette[i]=x;
        i++;
        TabPalette[i]=x;
        i++;
        TabPalette[i]=0;
        i++;
    }
    for (int x=0;x<1024;x++)
    {
        bfile.write((char*)&TabPalette[x], 1);
    }
}

/*
 * image given at the entry must be in the format : from bottom left to top
 */
void Bitmap::CreatePixel(ofstream &bfile, char * image)
{
    int white = 255;
    int black = 0;
    char pad = 0;
    int padToAdd = padWidth - width;
    for (int i=0; i<height; i++){
    	for(int j=0; j < width; j++)
            bfile.write((char*)(image+i*width+j), 1);
    	for(int j=0; j < padToAdd;j++)
    		bfile.write((char*)&pad,1);

    }
}
