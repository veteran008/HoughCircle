#pragma once
#include <string>
#include <vector>
//#include "../System/ViType.h"
#include <ViType.h>

#ifdef  DLL_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT __declspec(dllimport)   
#endif

enum ImageType
{
	BMP_TYPE = 0,
	TIFF_TYPE = 1,
	JPG_TYPE = 2,
	PNG_TYPE = 3,
	UNKNOWN_TYPE = 4
};

enum PIXEL_TYPE
{
	VIS_IMG_UINT8,
	VIS_IMG_UINT16,
	VIS_IMG_INT16,
	VIS_IMG_UINT32,
	VIS_IMG_INT32,
	VIS_IMG_REAL, //32
	VIS_IMG_DOUBLE, //64
};

enum IMAGE_TYPE
{
	VIS_IMG_GRAY = 1,
	VIS_IMG_RGB = 3,
	VIS_IMG_ARGB = 4
};

struct VIS_IMAGE_INFO
{
	IMG_INT width;
	IMG_INT height;
	IMG_INT byte_per_pixel;
	IMG_INT depth;
	IMG_INT page_num;
	IMG_INT channel;
	IMG_INT data_type;
	IMG_INT tiff_kind;
	IMG_INT planar_config;
	IMG_INT compression;
	IMG_INT unit_resolution;
	IMG_INT threshold;
	IMG_INT platte_size;
	IMG_INT pixel_size_x_mm;
	IMG_INT pixel_size_y_mm;
	IMG_INT fillorder;
	IMG_INT orientation;
	enum PIXEL_TYPE pixel_type;
	enum ImageType img_type;

	void clear()
	{
		memset(this, 0x00, sizeof(VIS_IMAGE_INFO));
	}

	VIS_IMAGE_INFO& operator=(const VIS_IMAGE_INFO& o) //÷ÿ–¥operator=
	{
		width = o.width;
		height = o.height;
		byte_per_pixel = o.byte_per_pixel;
		depth = o.depth;
		page_num = o.page_num;
		channel = o.channel;
		data_type = o.data_type;
		tiff_kind = o.tiff_kind;
		planar_config = o.planar_config;
		compression = o.compression;
		unit_resolution = o.unit_resolution;
		threshold = o.threshold;
		platte_size = o.platte_size;
		pixel_size_x_mm = o.pixel_size_x_mm;
		pixel_size_y_mm = o.pixel_size_y_mm;
		fillorder = o.fillorder;
		orientation = o.orientation;
		pixel_type = o.pixel_type;
		img_type = o.img_type;

		return *this;
	}

};

struct VIS_IMG_FILE_INFO
{
	IMG_UBYTE path[_MAX_PATH];
};

class  DLLEXPORT VisImage
{
public:
	VisImage();
	~VisImage();

	void GetImageInfo(VIS_IMAGE_INFO &imginfo);
	void SetImageInfo(VIS_IMAGE_INFO &imginfo);
	void SetImage(void **ppImage); //undone
	void GetImage(void **ppImage, IMG_INT channel=0);
	bool ReadImage(std::string filename, IMG_INT channel=0);
	bool WriteImage(std::string filename, IMG_INT channel=0); 
	void Reset();

private:
	IMG_INT Type(std::string filename);
	bool ExportColorPlane();
	IMG_INT ExportColorPlane(IMG_INT k, IMG_UBYTE **ppImage);
	IMG_INT ImportColorPlane(IMG_INT k, IMG_UBYTE **ppImage);

private:
	std::vector<IMG_UBYTE> m_vecPlatte;
	IMG_INT m_ImageType;
	VIS_IMAGE_INFO m_Imginfo;
	IMG_UBYTE *m_pCharData;
	std::string m_FileName;
	std::vector<void *> m_vecChannel;
};

