#pragma once
#include "ViType.h"

class VisBuf
{
public:
	VisBuf();
	~VisBuf();

	int set_IMG_UBBUF(IMG_UBBUF &srcBuf, IMG_UBYTE *ptr, IMG_SIZE size, IMG_UWORD  linestep);
	int set_IMG_WBUF(IMG_WBUF &srcBuf, IMG_WORD  *ptr, IMG_SIZE  size, IMG_UWORD linestep);
	int set_IMG_RBUF(IMG_RBUF &srcBuf, IMG_REAL  *ptr, IMG_SIZE  size, IMG_UWORD linestep);
private:

};