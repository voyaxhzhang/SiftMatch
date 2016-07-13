//Sift.h
/********************************************************************
	Sift
	created:	2015/04/23
	author:		LX_whu 
	purpose:	This file is for Sift function
	SIFTCPU:	copy from WuSift which  belongs to DPGRID 
	SIFTGPU:	copy from SIFTGPU which belongs to  Changchang Wu
*********************************************************************/
#if !defined Sift_h__LX_whu_2015_4_23
#define Sift_h__LX_whu_2015_4_23

typedef unsigned char	BYTE;
#define NORMALISE_IMAGE_SIDE(l) (l = l/2*2)	

#if !defined(ATMATCH_CODE) && !defined(SIFTMATCH_CODE)
	#ifdef SIFTMATCH_EXPORTS
		#define SIFTMATCH_API __declspec(dllexport)
	#else
		#define SIFTMATCH_API __declspec(dllimport)
		#if defined(_WIN64)  || defined(_X64)
			#pragma comment(lib,"SiftMatch_x64.lib") 
			#pragma message("Automatically linking with SiftMatch_x64.lib") 
		#else
			#ifdef _DEBUG_SIFTMATCH_API
				#pragma comment(lib,"SiftMatchD.lib") 
				#pragma message("Automatically linking with SiftMatchD.lib") 
			#else
				#pragma comment(lib,"SiftMatch.lib") 
				#pragma message("Automatically linking with SiftMatch.lib") 
			#endif
		#endif
	#endif
#else 
	#define SIFTMATCH_API
#endif

class  CSift
{
public:
	CSift();
	virtual ~CSift();
	SIFTMATCH_API virtual void	ParseParam(int argc, char **argv);
	SIFTMATCH_API bool			InitEnvi(bool bGPU = false);
	SIFTMATCH_API virtual int		RunSIFT(const BYTE* pImg, int nCols, int nRows);
	SIFTMATCH_API virtual	int		GetFeatureNum() const;
	/*
	****lcoation: x and y****
	****feature: scale and orientaion****
	****descriptors: 128****
	*/
	SIFTMATCH_API virtual void	GetFeatureVector(float* location, float* feature, float * descriptors) const;
	SIFTMATCH_API virtual void	GetFeatureVector(float* location, float* feature, BYTE * descriptors) const;
protected:
	virtual bool	_InitEnvi(){ return true; }
private:
	CSift*		m_sift;
};


typedef bool(*transform_err_fn)(const float* locL,const float* locR,const void* M,float dist);

class  CSiftMatch{
public:
	CSiftMatch();
	virtual ~CSiftMatch();
	SIFTMATCH_API bool			InitEnvi(bool bGPU = false);
	SIFTMATCH_API virtual void	SetMaxSift(int max_sift);
	SIFTMATCH_API virtual void	ParseParam(int argc, char **argv);

	//Specifiy descriptors to match, index = [0/1] for two features sets respectively
	//Option1, use float descriptors, and they be already normalized to 1.0
	SIFTMATCH_API virtual void SetDescriptors(int index, int num, const float* descriptors, int id = -1);
	//Option 2 unsigned char descriptors. They must be already normalized to 512
	SIFTMATCH_API virtual void SetDescriptors(int index, int num, const unsigned char * descriptors, int id = -1);
	//input feature location is a vector of [float x, float y, float skip[gap]]
	SIFTMATCH_API virtual void SetLocation(int index, const float* locations, int gap = 0);
	SIFTMATCH_API virtual void SetFeature(int index, const float* feature, int gap = 0);
public:
	/************************************************************************/
	/* all match_buffer below must  have the size >= max(ptNumL,ptNumR)	    */
	/************************************************************************/
	//match two sets of features, the function RETURNS the number of matches.
	//Given two normalized descriptor d1,d2, the distance here is acos(d1 *d2);
	SIFTMATCH_API virtual int  GetSiftMatch(
		int match_buffer[][2], //buffer to receive the matched feature indices
		float distmax = 0.7f,	//maximum distance of sift descriptor
		float ratiomax = 0.8f,	//maximum distance ratio
		int mutual_best_match = 1); //mutual best match or one way
	SIFTMATCH_API virtual int  GetNeighborSiftMatch(
		int match_buffer[][2], //buffer to receive the matched feature indices
		int maxNeighborNum = 2, //maximum number of neighbors of sift descriptor
		float distmax = 0.7f,	//maximum distance of sift descriptor
		int mutual_best_match = 1); 
	//use a guiding Homography H and a guiding Fundamental Matrix F to compute feature matches
	//the function returns the number of matches.
	SIFTMATCH_API virtual int  GetGuidedSiftMatch(
		int match_buffer[][2], //buffer to recieve 
		float H[3][3],			//homography matrix,  (Set NULL to skip)
		float F[3][3],			//fundamental matrix, (Set NULL to skip)
		float distmax = 0.7,	//maximum distance of sift descriptor
		float ratiomax = 0.8,   //maximum distance ratio
		float hdistmax = 32,    //threshold for |H * x1 - x2|_1 
		float fdistmax = 16,    //threshold for sampson error of x2'FX1
		int mutual_best_match = 1); //mutual best or one way
	SIFTMATCH_API virtual int  GetXformSiftMatch(
		int match_buffer[][2],
		transform_err_fn er_fn, void* xform, float errdistmax,
		float distmax = 0.7,
		float ratiomax = 0.8,		
		int mutual_best_match = 1
		);
	//match two sets of feature, and the matches is limited to scale/orientation
	//if ratiomin_scale<0 || ratiomax_scale < 0, skip scale feature
	//if diffmin_ori<0 || diffmax_ori > PI ,skip orientation feature
	//diff_scale = scale_right/scale_left
	//diff_orietation = |orietation_right-orietation_left|
	SIFTMATCH_API virtual	int		GetLimitedSiftMatch(
		int match_buffer[][2],
		float ratiomin_scale,float ratiomax_scale,
		float diffmin_ori,float diffmax_ori,
		float distmax = 0.7f,
		float ratiomax = 0.8f,
		int mutual_best_match = 1
		);
	SIFTMATCH_API virtual int		GetLimitedGuidedSiftMatch(
		int match_buffer[][2],
		float H[3][3],
		float F[3][3],
		float ratiomin_scale, float ratiomax_scale,
		float diffmin_ori, float diffmax_ori,
		float distmax = 0.7f,
		float ratiomax = 0.8f,
		float hdistmax = 32,
		float fdistmax = 16,
		int mutual_best_match = 1
		);
	SIFTMATCH_API virtual int		GetLimitedXformSiftMatch(
		int match_buffer[][2],
		transform_err_fn er_fn, void* xform, float errdistmax,
		float ratiomin_scale, float ratiomax_scale,
		float diffmin_ori, float diffmax_ori,		
		float distmax = 0.7f,
		float ratiomax = 0.8f,
		int mutual_best_match = 1
		);	
protected:
	virtual bool	_InitEnvi(){ return true; }
protected:
	int				m_siftmatch_max_num;
private:
	CSiftMatch*		m_sift_matcher;
	int				m_match_condition[2];
};

SIFTMATCH_API bool		CheckLimitedScale(float ratiomin_scale, float ratiomax_scale);
SIFTMATCH_API bool		CheckLimitedOrientation(float diffmin_ori, float diffmax_ori);

//decrease search radius step by step. 
//three step: 1)GetSiftMatch 2)calculate H and F 3)GetGuidedSiftMatch
SIFTMATCH_API int		GetComboSiftMatch(
	const float* locationsL, const BYTE* desL, int ptNumL,
	const float* locationsR, const BYTE* desR, int ptNumR,
	int match_buffer[][2],
	bool bUseH = true, bool bUseF = true,
	float distmax = 0.7f,
	float ratiomax = 0.8f,
	int mutual_best_match = 1,
	bool bGPU = false);
SIFTMATCH_API int		GetSimilaritySiftMatch(
	const float* locationsL, const BYTE* desL, const float* featureL, int ptNumL,
	const float* locationsR, const BYTE* desR, const float* featureR, int ptNumR,
	int match_buffer[][2],
	float distmax = 0.7f,
	float ratiomax = 0.8f,
	int mutual_best_match = 1,
	bool bGPU = false);

#ifndef _FPT4D
#define _FPT4D
typedef struct tagFPT4D
{
	float xl, yl, xr, yr;
}FPT4D;
#endif

SIFTMATCH_API FPT4D*	 GetSiftMatch(
	int* ptSum,
	const BYTE *pImgL, int colsL, int rowsL,
	const BYTE *pImgR, int colsR, int rowsR,
	bool bUseH = false,bool bUseF = false,
	float distmax = 0.7f,
	float ratiomax = 0.8f,
	int mutual_best_match = 1,
	bool bGPU = false
	);
SIFTMATCH_API FPT4D*	 GetNeighborSiftMatch(
	int* ptSum,
	const BYTE *pImgL, int colsL, int rowsL,
	const BYTE *pImgR, int colsR, int rowsR,
	int maxNeighborNum = 2, 
	float distmax = 0.5f,	
	int mutual_best_match = 1
	);

SIFTMATCH_API FPT4D*	 GetXformSiftMatch(
	int* ptSum,
	const BYTE *pImgL, int colsL, int rowsL,
	const BYTE *pImgR, int colsR, int rowsR,
	transform_err_fn er_fn, void* xform, float errdistmax,
	float distmax = 0.7f,
	float ratiomax = 0.8f,
	int mutual_best_match = 1
	);

SIFTMATCH_API int		ExtractSiftFromBuffer(const BYTE* pImg, int nCols, int nRows, float** loc, float** fea, BYTE** des, bool bGPU);
SIFTMATCH_API void		FreeSiftPtsMem(void* pt);
#endif // Sift_h__LX_whu_2015_4_23