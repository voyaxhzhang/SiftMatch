#include "SiftMatch.h"

#include "algorithms/Transformation.h"

#ifdef _FLANN_CPP
#include "flann/flann.hpp"
#else
#include "flann/flann.h"
#pragma comment(lib,"flann.lib") 
#endif

#ifdef _SIFTGPU

#ifdef _WIN32
#include <windows.h>
#define FREE_MYLIB FreeLibrary
#define GET_MYPROC GetProcAddress
#else
#include <dlfcn.h>
#define FREE_MYLIB dlclose
#define GET_MYPROC dlsym
#endif

#include "src/SiftGPU/SiftGPU.h"
#include "include/GL/glew.h"

#if !defined(SIFTGPU_STATIC) && !defined(SIFTGPU_DLL_RUNTIME) 
// SIFTGPU_STATIC comes from compiler
#define SIFTGPU_DLL_RUNTIME
// Load at runtime if the above macro defined
// comment the macro above to use static linking
#endif

#ifndef SIFTGPU_DLL_RUNTIME

#define SIFTGPU_DLL
#if defined(_WIN64)  || defined(_X64)
#pragma comment(lib, "siftgpu64.lib")
#else
#ifdef _DEBUG 
#pragma comment(lib, "siftgpu_d.lib")
#else
#pragma comment(lib, "siftgpu.lib")
#endif
#endif

#endif

#endif
#include "../MatchBase/GlobalUtil_AtMatch.h"

#ifdef _SIFTGPU
struct SiftGPU_dll_register
{
	SiftGPU_dll_register()
	{
#ifdef SIFTGPU_DLL_RUNTIME
		hsiftgpu = NULL;		
#endif
		b_load_fun = false;
		pCreateNewSiftGPU = NULL;
		pCreateNewSiftMatchGPU = NULL;
	}
	~SiftGPU_dll_register()
	{
#ifdef SIFTGPU_DLL_RUNTIME
		if (hsiftgpu) FREE_MYLIB(hsiftgpu);
#endif
	}
	bool initial_dll()
	{
		if (b_load_fun) return true;
#ifdef SIFTGPU_DLL_RUNTIME
	#ifdef _WIN32
#if defined(_WIN64)  || defined(_X64)
		HMODULE  hsiftgpu = LoadLibrary("siftgpu64.dll");
#else
		#ifdef _DEBUG
				HMODULE  hsiftgpu = LoadLibrary("siftgpu_d.dll");
		#else
				HMODULE  hsiftgpu = LoadLibrary("siftgpu.dll");
		#endif
#endif
	#else
			void * hsiftgpu = dlopen("libsiftgpu.so", RTLD_LAZY);
	#endif
			if (hsiftgpu == NULL) {
				GlobalParam::g_sift_extract_gpu = false;
				GlobalParam::g_sift_match_gpu = false;
				return false;
			}
		pCreateNewSiftGPU = (SiftGPU* (*) (int)) GET_MYPROC(hsiftgpu, "CreateNewSiftGPU");
		pCreateNewSiftMatchGPU = (SiftMatchGPU* (*)(int)) GET_MYPROC(hsiftgpu, "CreateNewSiftMatchGPU");
#else
		pCreateNewSiftGPU = CreateNewSiftGPU;
		pCreateNewSiftMatchGPU = CreateNewSiftMatchGPU;
#endif
		b_load_fun = true;
		return true;
	}
	SiftGPU* CreateSiftGPU(int np = 1){
		if (!initial_dll()) return NULL;
		return pCreateNewSiftGPU(np);
	}
	SiftMatchGPU* CreateSiftMatchGPU(int max_sift = 4096)
	{
		if (!initial_dll()) return NULL;
		return pCreateNewSiftMatchGPU(max_sift);
	}
#ifdef SIFTGPU_DLL_RUNTIME	
	#ifdef _WIN32
		HMODULE hsiftgpu;
	#else
		void* hsiftgpu;
	#endif
#endif
	bool	b_load_fun;
	SiftGPU* (*pCreateNewSiftGPU)(int);
	SiftMatchGPU* (*pCreateNewSiftMatchGPU)(int);
} g_siftgpu_dll_register;

#endif // _SIFTGPU

#ifdef WUSIFT_V1
#include "DPGRID/WuSift-v1.0.h"
#else
#include "DPGRID/WuSift-v2.0.h"
#endif

#include <vector>
using namespace std;

#define PI 3.1415926

#ifdef _SIFTGPU
class CSiftGPU : public CSift
{
public:
	CSiftGPU();
	virtual ~CSiftGPU();
	virtual void	ParseParam(int argc, char **argv);
	virtual int		RunSIFT(const BYTE* pImg, int nCols, int nRows);
	virtual	int		GetFeatureNum() const;
	virtual void	GetFeatureVector(float* location, float* feature, float * descriptors) const;
	virtual void	GetFeatureVector(float* location, float* feature, BYTE * descriptors) const;
public:
	virtual bool	_InitEnvi();
private:
	SiftGPU*	m_sift;
};
#endif

class CSiftCPU : public CSift,public CWuSift
{
public:
	CSiftCPU();
	virtual ~CSiftCPU();
	virtual void	ParseParam(int argc, char **argv);
	
	virtual int		RunSIFT(const BYTE* pImg, int nCols, int nRows);
	virtual	int		GetFeatureNum() const;
	virtual void	GetFeatureVector(float* location, float* feature, float * descriptors) const;
	virtual void	GetFeatureVector(float* location, float* feature, BYTE * descriptors) const;
public:
	virtual bool	_InitEnvi();
private:
	int				m_fea_num;
	FEATDES*		m_fea_buf;
	int				m_zoom_ratio;
};

CSift::CSift(){
	m_sift = NULL;
}

CSift::~CSift(){
	if (m_sift) delete m_sift;
}

void* CSift::operator new (size_t  size){
	void * p = malloc(size);
	if (p == 0)
	{
		const std::bad_alloc ba;
		throw ba;
	}
	return p;
}

bool CSift::InitEnvi(bool bGPU)
{
	if (m_sift)  {
		delete m_sift; m_sift = NULL;
	}
#ifdef _SIFTGPU
	if (bGPU) {
		m_sift = new CSiftGPU;
		if (m_sift){
			if (m_sift->_InitEnvi()){
				LogPrint(0, "[USE SIFTGPU]");
				return true;
			}
			delete m_sift;	m_sift = NULL;
		}
		LogPrint(0, "SIFTGPU is not supported.Check if exist SIFTGPU library.");
	}
#endif
	LogPrint(0, "[USE SIFTCPU]");
	m_sift = new CSiftCPU;
	return m_sift ? m_sift->_InitEnvi() : false;
}

void CSift::ParseParam(int argc, char **argv)
{
	if (!m_sift&&!InitEnvi(GlobalParam::g_sift_extract_gpu)) return;
	return m_sift->ParseParam(argc, argv);
}

int	 CSift::RunSIFT(const BYTE* pImg, int nCols, int nRows)
{
	if (!m_sift&&!InitEnvi(GlobalParam::g_sift_extract_gpu)) return 0;
	return m_sift->RunSIFT(pImg, nCols, nRows);
}

int CSift::GetFeatureNum() const{
	if (!m_sift) return 0;
	return m_sift->GetFeatureNum();
}

void CSift::GetFeatureVector(float* location, float* feature, float * descriptors) const
{
	if (!m_sift) return;
	m_sift->GetFeatureVector(location, feature, descriptors);
}

void CSift::GetFeatureVector(float* location, float* feature, BYTE * descriptors) const
{
	if (!m_sift) return;
	m_sift->GetFeatureVector(location, feature, descriptors);
}

#ifdef _SIFTGPU
CSiftGPU::CSiftGPU(){
	m_sift = g_siftgpu_dll_register.CreateSiftGPU(1);
}

CSiftGPU::~CSiftGPU(){
	if (m_sift) {
		m_sift->DestroyContextGL();
		delete m_sift;
	}
}

bool CSiftGPU::_InitEnvi()
{
	if (!m_sift) return false;
	if (!GlobalParam::g_siftgpu_initial) {
		if (m_sift->CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED){
			LogPrint(ERR_ATMCH_FLAG, "Fail to initialize SIFTGPU!");
			return false;
		}
		GlobalParam::g_siftgpu_initial = true;
	}
	else
		m_sift->VerifyContextGL();
//	GlobalParam::g_max_buffer_size = 4 * 1024 * 1024;
	return true;
}

void CSiftGPU::ParseParam(int argc, char **argv)
{
	m_sift->ParseParam(argc, argv);
}

int	 CSiftGPU::RunSIFT(const BYTE* pImg, int nCols, int nRows)
{
	return m_sift->RunSIFT(nCols, nRows, pImg, GL_LUMINANCE, GL_UNSIGNED_BYTE);
}

int CSiftGPU::GetFeatureNum() const
{
	return m_sift->GetFeatureNum();
}


void CSiftGPU::GetFeatureVector(float* location, float* feature, float * descriptors) const
{
	int num = GetFeatureNum();
	if (num < 1) return;
	vector<SiftGPU::SiftKeypoint> keys;		keys.resize(num);
	m_sift->GetFeatureVector(&keys[0], descriptors);
	for (int i = 0; i < num; i++){
		float* pBuf = (float*)&keys[i];
		memcpy(location, pBuf, sizeof(float)* 2);
		memcpy(feature, pBuf + 2, sizeof(float)* 2);
		location += 2;
		feature += 2;
	}
}

void CSiftGPU::GetFeatureVector(float* location, float* feature, BYTE * descriptors) const
{
	int num = GetFeatureNum();
	if (num < 1) return;
	vector<float> des;	des.resize(128 * num);
	float* pDes = &des[0];
	GetFeatureVector(location, feature, pDes);
	
	for (int i = 0; i < num; i++){
		CvtDescriptorsf2uc(pDes, descriptors, 128);
		pDes += 128;
		descriptors += 128;
	}
}
#endif

CSiftCPU::CSiftCPU(){
	m_fea_buf = NULL;
	m_zoom_ratio = 1;
}

CSiftCPU::~CSiftCPU(){
	if (m_fea_buf) Free_Feat(m_fea_buf);
}

bool CSiftCPU::_InitEnvi()
{
	return true;
}

void CSiftCPU::ParseParam(int argc, char **argv)
{

}

//parallel calculate
int CSiftCPU::RunSIFT(const BYTE* pImg, int nCols, int nRows)
{
	if (m_fea_buf) Free_Feat(m_fea_buf);
	int miSz = nCols < nRows ? nCols : nRows;
	if (miSz>1024)  miSz = miSz>>1;
	m_fea_num = 0;
	m_fea_buf = (FEATDES*)(Extract_Feat(&m_fea_num, (BYTE*)pImg, nCols, nRows, &m_zoom_ratio, true, miSz));
	return m_fea_num;
}

int CSiftCPU::GetFeatureNum() const{
	return m_fea_num;
}

#ifdef WUSIFT_V1
void CSiftCPU::GetFeatureVector(float* location, float* feature, float * descriptors) const
{
	if (!m_fea_buf) return;
	FEATDES* pFea = m_fea_buf;
	for (int i = 0; i < m_fea_num; i++, pFea++)
	{
		*location++ = (float)(pFea->x*m_zoom_ratio);
		*location++ = (float)(pFea->y*m_zoom_ratio);
		*feature++ = pFea->scl*m_zoom_ratio;
		*feature++ = pFea->ori;
		
		int v;
		float *pd1,*pd2, sq = 0;
		//int v;
		//normalize
		pd1 = pFea->descr;
		for (v = 0; v < 128; v++, pd1++)	sq += (*pd1)*(*pd1);
		sq = 1.0f / sqrtf(sq);
		//truncate to .2
		pd1 = pFea->descr;	pd2 = descriptors;
		for (v = 0; v < 128; v++, pd1++, pd2++)	*pd2 = min(*pd1*sq, 0.2f);

		//renormalize
		pd2 = descriptors; sq = 0;
		for (v = 0; v < 128; v++, pd2++)	sq += (*pd2)*(*pd2);
		sq = 1.0f / sqrtf(sq);

		for (v = 0; v < 128; v++, descriptors++)	*descriptors = *descriptors*sq;

	}
}

void	CSiftCPU::GetFeatureVector(float* location, float* feature, BYTE * descriptors) const{
	if (!m_fea_buf || m_fea_num<1) return;
	vector<float> des;	des.resize(128 * m_fea_num);
	float* pDes = &des[0];
	GetFeatureVector(location, feature, pDes);

	for (int i = 0; i < m_fea_num; i++){
		CvtDescriptorsf2uc(pDes, descriptors, 128);
		pDes += 128;
		descriptors += 128;
	}
}
#else
void CSiftCPU::GetFeatureVector(float* location, float* feature, float * descriptors) const
{
	if (!m_fea_buf) return;
	FEATDES* pFea = m_fea_buf;
	for (int i = 0; i < m_fea_num; i++, pFea++)
	{
		*location++ = (float)(pFea->x*m_zoom_ratio);
		*location++ = (float)(pFea->y*m_zoom_ratio);
		*feature++ = pFea->scl*m_zoom_ratio;
		*feature++ = pFea->ori;
		BYTE* pd = pFea->descr;
		for (int j = 0; j < 128; j++ ){
			*descriptors++ = *pd++ / 512.0f;
		}
	}
}
void	CSiftCPU::GetFeatureVector(float* location, float* feature, BYTE * descriptors) const{
	if (!m_fea_buf) return;
	FEATDES* pFea = m_fea_buf;
	for (int i = 0; i < m_fea_num; i++, pFea++)
	{
		*location++ = (float)(pFea->x*m_zoom_ratio);
		*location++ = (float)(pFea->y*m_zoom_ratio);
		*feature++ = pFea->scl*m_zoom_ratio;
		*feature++ = pFea->ori;
		memcpy(descriptors, pFea->descr, sizeof(BYTE)* 128);		
		descriptors += 128;
	}
}
#endif

typedef bool(*check_feature_fn)(float sl, float sr, float minS,  float maxS);

bool check_scale(float sl,  float sr,  float minS, float maxS){
	float s = sr / sl;
	if (minS <= s && s <= maxS) return true;
	return false;
}
bool check_orientation(float sl, float sr, float minS, float maxS){
	float s = fabs(sr - sl);	if (s>PI) s = 2 * PI - s;
	if (minS < s && s < maxS) return true;
	return false;
}

bool homog_check_err(const float* feaL, const  float* feaR, const  void* M, float dist){
	float to[2];
	GeoTransform::homog_xform_pt(feaL, to, (const float*)M);
	float er[2];
	er[0] = fabs(to[0] - feaR[0]);
	er[1] = fabs(to[1] - feaR[1]);
	if (dist>er[0] && dist>er[1]) return true;
	return false;
}

bool fundamental_check_err(const float* feaL, const float* feaR, const  void* M, float dist){
	if (GeoTransform::fundamental_xfer_distance(feaL, feaR, (const float*)M) < dist) return true;
	return false;
}

/************************************************************************/
/* flann use random kdtree, so result may change in different run time  */
/************************************************************************/
class CSiftMatchCPU : public CSiftMatch
{
private:
#ifdef _FLANN_CPP
	typedef flann::L2<BYTE>		Des_Distance;
	typedef flann::L2<float>	Loc_Distance;
	typedef flann::Index<Des_Distance>	Des_Index;
	typedef flann::Index<Loc_Distance>	Loc_Index;
	Des_Index*				m_flann_descriptors_index[2];
	Loc_Index*				m_flann_location_index[2];
#else
	struct FLANNParameters		m_flann_parameters;
	flann_index_t				m_flann_descriptors_index[2];
	flann_index_t				m_flann_location_index[2];
#endif
	int						m_sift_num[2];
	vector<BYTE>			m_descriptors[2];
	vector<float>			m_location[2];
	vector<float>			m_feature[2];
	int						m_flann_descriptors_checks;
	int						m_flann_descriptors_trees;
protected:
	inline bool	BuildDescriptorsIndex(int index){
		LogPrint(0, "Build flann index for [%d]descriptors...", index);
#ifdef _FLANN_CPP
		if (m_flann_descriptors_index[index]) delete m_flann_descriptors_index[index];	
		m_flann_descriptors_index[index] = new Des_Index(flann::Matrix<BYTE>(m_descriptors[index].data(), m_descriptors[index].size() / 128, 128),flann::KDTreeIndexParams(m_flann_descriptors_trees) );//
		m_flann_descriptors_index[index]->buildIndex();
#else
		if (m_flann_descriptors_index[index]) flann_free_index(m_flann_descriptors_index[index], &m_flann_parameters);
		float speedup;
		m_flann_descriptors_index[index] = flann_build_index_byte(m_descriptors[index].data(), m_descriptors[index].size() / 128, 128, &speedup, &m_flann_parameters);
#endif
		return true;
	}
	inline bool	BuildLocationIndex(int index){
		LogPrint(0, "Build flann index for [%d]location...", index);
#ifdef _FLANN_CPP
		if (m_flann_location_index[index])  delete m_flann_location_index[index];
		m_flann_location_index[index] = new Loc_Index(flann::Matrix<float>(m_location[index].data(), m_location[index].size() / 2, 2), flann::KDTreeIndexParams(4));
		m_flann_location_index[index]->buildIndex();
#else
		if (m_flann_location_index[index]) flann_free_index(m_flann_location_index[index], &m_flann_parameters);
		float speedup;
		m_flann_location_index[index] = flann_build_index_float(m_location[index].data(), m_location[index].size() / 2, 2, &speedup, &m_flann_parameters);
#endif
		return true;
	}
	inline int FindDesNearestNeighborsIndex(void* index_id,
		unsigned char* testset,
		int trows,
		int* indices,
		float* dists,
		int nn
		){
#ifdef  _FLANN_CPP
		return ((Des_Index*)index_id)->knnSearch(flann::Matrix<BYTE>(testset, trows, 128), flann::Matrix<int>(indices, trows, nn), flann::Matrix<float>(dists, trows, nn), nn,flann::SearchParams(m_flann_descriptors_checks) );//
#else
		return flann_find_nearest_neighbors_index_byte(index_id, testset, trows, indices, dists, nn, &m_flann_parameters);
#endif
	}
	inline int RadiusSearch(
		void* index_ptr, 
		unsigned char* query, 
		int* indices,
		float* dists, 
		int max_nn, 
		float radius
	){
#ifdef  _FLANN_CPP
		flann::SearchParams params;	
		params.checks = m_flann_descriptors_checks;
		params.sorted = 1;
		params.max_neighbors = max_nn;
		return ((Des_Index*)index_ptr)->radiusSearch(flann::Matrix<BYTE>(query, 1, 128), flann::Matrix<int>(indices, 1, max_nn), flann::Matrix<float>(dists, 1, max_nn), radius, params);
#else
	return flann_radius_search_byte(index_ptr, query, indices, dists, max_nn, radius, &m_flann_parameters);
#endif
	}
public:
	CSiftMatchCPU(){
		m_sift_num[0] = m_sift_num[1] = 0;		

		m_flann_descriptors_index[0] = m_flann_descriptors_index[1] = NULL;
		m_flann_location_index[0] = m_flann_location_index[1] = NULL;
		m_flann_descriptors_trees = 8;
		m_flann_descriptors_checks = 128;
#ifndef _FLANN_CPP		
		m_flann_parameters = DEFAULT_FLANN_PARAMETERS;
		m_flann_parameters.algorithm = FLANN_INDEX_KDTREE;
		m_flann_parameters.trees = m_flann_descriptors_trees;
		m_flann_parameters.log_level = FLANN_LOG_INFO;
		m_flann_parameters.checks = m_flann_descriptors_checks;
		m_flann_parameters.sorted = 1;
#endif
	}
	virtual ~CSiftMatchCPU(){
#ifdef _FLANN_CPP
		if (m_flann_descriptors_index[0]) delete m_flann_descriptors_index[0];
		if (m_flann_descriptors_index[1]) delete m_flann_descriptors_index[1];
		if (m_flann_location_index[0]) delete m_flann_location_index[0];
		if (m_flann_location_index[1]) delete m_flann_location_index[1];
#else 
		if (m_flann_descriptors_index[0]) flann_free_index(m_flann_descriptors_index[0], &m_flann_parameters);
		if (m_flann_descriptors_index[1]) flann_free_index(m_flann_descriptors_index[1], &m_flann_parameters);
		if (m_flann_location_index[0]) flann_free_index(m_flann_location_index[0], &m_flann_parameters); 
		if (m_flann_location_index[1]) flann_free_index(m_flann_location_index[1], &m_flann_parameters); 
#endif
	}
	virtual void	SetMaxSift(int max_sift){
		m_siftmatch_max_num = max_sift;
	}
	virtual bool	_InitEnvi()
	{
		return true;
	}
	virtual void	ParseParam(int argc, char **argv)
	{
		for (int i=0; i<argc; i++ ){
			if (!strcmp(argv[i],"-trees")){
				if (++i >= argc) break;
				int t = atoi(argv[i]);
				if (t >= 2 && t <= 128) m_flann_descriptors_trees = (BYTE)t;
			}
			else if (!strcmp(argv[i], "-checks")){
				if (++i >= argc) break;
				int t = atoi(argv[i]);
				if (t >= 2 && t <= 128) m_flann_descriptors_checks = (BYTE)t;
			}
		}
#ifndef _FLANN_CPP
		m_flann_parameters.checks = m_flann_descriptors_checks;
		m_flann_parameters.trees = m_flann_descriptors_trees;
#endif
	}
	virtual void	SetDescriptors(int index, int num, const unsigned char * descriptors, int id /* = -1 */)
	{
		m_descriptors[index].resize(num*128);
		memcpy(m_descriptors[index].data(), descriptors, num*sizeof(BYTE)* 128);
		BuildDescriptorsIndex(index);
		m_sift_num[index] = num;
	}
	virtual void	SetDescriptors(int index, int num, const float * descriptors, int id){
		
		m_descriptors[index].resize(num * 128);
		
		BYTE*			pDst = &m_descriptors[index][0];
		const float*	pSrc = descriptors;

		for (int i = 0; i < num; i++){
			CvtDescriptorsf2uc(pSrc, pDst, 128);
			pSrc += 128;	pDst += 128;
		}

		BuildDescriptorsIndex(index);
		m_sift_num[index] = num;
	}
	virtual void	SetLocation(int index, const float* locations, int gap)
	{
		int num = m_sift_num[index];
		m_location[index].resize(num*2);
		const float* pSrc = locations;
		float* pDst = &m_location[index][0];
		for (int i = 0; i < num; i++){
			memcpy(pDst, pSrc, 2 * sizeof(float));
			pDst += 2;	pSrc += 2 + gap;
		}
//		BuildLocationIndex(index);
	}
	virtual void	SetFeature(int index, const float* feature, int gap)
	{
		int num = m_sift_num[index];
	
		m_feature[index].resize(num*2);
		const float* pSrc = feature;
		float* pDst = m_feature[index].data();
		for (int i = 0; i < num; i++){
			memcpy(pDst, pSrc, 2 * sizeof(float));
			pDst += 2;
			pSrc += 2 + gap;
		}
	}
	virtual int		GetSiftMatch(int match_buffer[][2], float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
	{
		Prepare(&distmax, &ratiomax);

		int* r_in_l = new int[m_sift_num[1]];
		LogPrint(0, "Find nearest neighbors ...");
		int num_rinl = GetBestMatch(r_in_l, m_flann_descriptors_index[0], &m_descriptors[1][0], m_sift_num[1], distmax, ratiomax);
		LogPrint(0, "Get correspond num = %d", num_rinl);
		if (!mutual_best_match){
			int* pR2L = r_in_l;
			for (int i = 0; i < m_sift_num[1]; i++,pR2L++ ){
				if (*pR2L < 0) continue;
				(*match_buffer)[0] = *pR2L;
				(*match_buffer)[1] = i;
				match_buffer++;
			}
			delete r_in_l;
			return num_rinl;
		}
		LogPrint(0, "Find nearest neighbors in negative direction...");
		int* l_in_r = new int[m_sift_num[0]];
		int num_linr = GetBestMatch(l_in_r, m_flann_descriptors_index[1], &m_descriptors[0][0], m_sift_num[0], distmax, ratiomax);
		LogPrint(0, "Get correspond num = %d", num_linr);
		int* pR2L = r_in_l;	int cnt = 0;

		for (int i = 0; i < m_sift_num[1]; i++, pR2L++){
			if (*pR2L < 0) continue;
			if (l_in_r[*pR2L] != i) continue;
			(*match_buffer)[0] = *pR2L;
			(*match_buffer)[1] = i;
			match_buffer++;
			cnt++;
		}
		delete r_in_l;
		delete l_in_r;

		LogPrint(0, "final correspond num = %d", cnt);
		return cnt;
	}
	virtual int		GetNeighborSiftMatch(int match_buffer[][2], int maxNeighborNum /* = 2 */, float distmax /* = 0.7f */, int mutual_best_match /* = 1 */);
	virtual int		GetGuidedSiftMatch(int match_buffer[][2], float H[3][3], float F[3][3], float distmax, float ratiomax, float hdistmax, float fdistmax, int mutual_best_match)
	{
		Prepare(&distmax, &ratiomax);

		
		transform_err_fn er_fn[2] = { homog_check_err, fundamental_check_err };
		float	trans_dist[2] = { hdistmax, fdistmax };
		float*	M[2] = { H[0], F[0] };
		return	GetConditionMatch(match_buffer, distmax, ratiomax, er_fn, (void**)M, trans_dist, 2, NULL, 0, 0, NULL, 0, 0, mutual_best_match);
	}
	virtual int GetXformSiftMatch(
		int match_buffer[][2],
		transform_err_fn er_fn, void* xform, float errdistmax,
		float distmax ,
		float ratiomax,
		int mutual_best_match
		)
	{
		Prepare(&distmax, &ratiomax);
		return	GetConditionMatch(match_buffer, distmax, ratiomax, &er_fn, &xform, &errdistmax, 1, NULL, 0, 0, NULL, 0, 0, mutual_best_match);
	}
	virtual int		GetLimitedSiftMatch(int match_buffer[][2], float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
	{
		Prepare(&distmax, &ratiomax);
		return GetConditionMatch(match_buffer, distmax, ratiomax, NULL, NULL, NULL, 0,
			check_scale, ratiomin_scale, ratiomax_scale,
			check_orientation, diffmin_ori, diffmax_ori,
			mutual_best_match);
	}
	virtual int		GetLimitedGuidedSiftMatch(int match_buffer[][2], float H[3][3], float F[3][3], float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax, float ratiomax, float hdistmax, float fdistmax, int mutual_best_match)
	{
		Prepare(&distmax, &ratiomax);
		transform_err_fn er_fn[2] = { homog_check_err, fundamental_check_err };
		float	trans_dist[2] = { hdistmax, fdistmax };
		float*	M[2] = { H[0], F[0] };
		return	GetConditionMatch(match_buffer, distmax, ratiomax,er_fn,(void**)M,trans_dist,2,
			check_scale, ratiomin_scale, ratiomax_scale,
			check_orientation, diffmin_ori, diffmax_ori,
			mutual_best_match);
	}
	virtual int		GetLimitedXformSiftMatch(int match_buffer[][2], transform_err_fn er_fn, void* xform, float errdistmax, float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax, float ratiomax, int mutual_best_match)
	{
		Prepare(&distmax, &ratiomax);
		return	GetConditionMatch(match_buffer, distmax, ratiomax, &er_fn, &xform, &errdistmax, 1,
			check_scale, ratiomin_scale, ratiomax_scale,
			check_orientation, diffmin_ori, diffmax_ori,
			mutual_best_match);
	}
protected:
	void		Prepare(float* distmax , float* ratiomax ){
		if (distmax){
			*distmax = 2 * sin(*distmax / 2) * 512;	*distmax = *distmax**distmax;
		}
		if(ratiomax) *ratiomax = *ratiomax**ratiomax;
	}
	int		GetConditionMatch(
		int match_buffer[][2], 
		float desdistmax, float ratiomax,
		transform_err_fn* er_fn, void** M, float* distmax, int num_fn,
		check_feature_fn chek_scale_fn, float minS, float maxS,
		check_feature_fn chek_ori_fn, float minO, float maxO,
		int mutual_best_match
		){
		int i;
		int* r_in_l = new int[m_sift_num[1]];
		LogPrint(0, "Find guided nearest neighbors ...");
		BYTE*	pDes1;
		BYTE*	pDes2 = &m_descriptors[1][0];
		float*	pLoc2 = &m_location[1][0];
		float*	pLoc1 = &m_location[0][0];

		float*	pScl1 = NULL;
		float*  pScl2 = NULL;
		float*	pOri1 = NULL;
		float*	pOri2 = NULL;
		int		scl_buf_space = 0;
		int		ori_buf_space = 0;

		if (chek_scale_fn){
			pScl1 = &m_feature[0][0];
			pScl2 = &m_feature[1][0];
			scl_buf_space = 2;
		}
		if (chek_ori_fn){
			pOri1 = &m_feature[0][1];
			pOri2 = &m_feature[1][1];
			ori_buf_space = 2;
		}


		int num_rinl;
		num_rinl = GetBestConditionMatch(r_in_l, pDes2, pLoc2, pScl2, pOri2, m_sift_num[1],
			m_flann_descriptors_index[0], desdistmax, ratiomax,
			er_fn, pLoc1, 2, M, distmax, num_fn, false,
			chek_scale_fn, pScl1, scl_buf_space, minS, maxS,
			chek_ori_fn, pOri1, ori_buf_space, minO, maxO);

// 		for (i = 0; i < m_sift_num[1]; i++){
// 			r_in_l[i] = GetPointMatch(pDes2, pLoc2, pScl2, pOri2, m_flann_descriptors_index[0], desdistmax, ratiomax,
// 				er_fn, pLoc1, 2, M, distmax, num_fn, false,
// 				chek_scale_fn, pScl1, 2, minS, maxS,
// 				chek_ori_fn, pOri1, 2, minO, maxO);
// 			pDes2 += 128;
// 			pLoc2 += 2;
// 			pScl2 += scl_buf_space;
// 			pOri2 += ori_buf_space;
// 			if (r_in_l[i]>0){
// 				match_buffer[num_rinl][0] = r_in_l[i];
// 				match_buffer[num_rinl][1] = i;
// 				num_rinl++;
// 			}
// 		}

		LogPrint(0, "Get correspond num = %d", num_rinl);
		if (!mutual_best_match || num_rinl < 1){
			if (num_rinl > 0){
				num_rinl = 0;
				for (i = 0; i < m_sift_num[1]; i++){
					if (r_in_l[i]>=0){
						match_buffer[num_rinl][0] = r_in_l[i];
						match_buffer[num_rinl][1] = i;
						num_rinl++;
					}
				}
			}
			delete r_in_l;
			
			return num_rinl;
		}

		pDes1 = new BYTE[num_rinl * 128];
		pLoc1 = new float[num_rinl * 2];
		if (pScl1 || pOri1)
		{
			pScl1 = new float[num_rinl * 2];
			pOri1 = pScl1 + 1;
		}
		
		num_rinl = 0;
		for (i = 0; i < m_sift_num[1]; i++){
			if (r_in_l[i]>=0){
				match_buffer[num_rinl][0] = r_in_l[i];
				match_buffer[num_rinl][1] = i;
				memcpy(pDes1 + 128 * num_rinl, &m_descriptors[0][128 * r_in_l[i]], 128 );
				memcpy(pLoc1 + 2 * num_rinl, &m_location[0][2 * r_in_l[i]], 2 * sizeof(float));
				if (pScl1) pScl1[2*num_rinl] = m_feature[0][2 * r_in_l[i]];
				if (pOri1) pOri1[2*num_rinl] = m_feature[0][2 * r_in_l[i] + 1];

				num_rinl++;
			}
		}
		delete r_in_l;

		pLoc2 = &m_location[1][0];
		if (chek_scale_fn){
			pScl2 = &m_feature[1][0];
		}
		if (chek_ori_fn){
			pOri2 = &m_feature[1][1];
		}

// 		int mn = 0, idx;
// 		for (i = 0; i < num_rinl; i++){
// 			int idx1 = match_buffer[i][0];
// 			pDes1 = &m_descriptors[0][0] + idx1 * 128;
// 			pLoc1 = &m_location[0][0] + idx1 * 2;
// 			if (chek_scale_fn) pScl1 = &m_feature[0][0] + idx1 * 2;
// 			if (chek_ori_fn)	pOri1 = &m_feature[0][1] + idx1 * 2;
// 			idx = GetPointMatch(pDes1, pLoc1, pScl1, pOri1, m_flann_descriptors_index[1], desdistmax, ratiomax,
// 				er_fn, pLoc2, 2, M, distmax, num_fn,true,
// 				chek_scale_fn, pScl2, 2, minS, maxS,
// 				chek_ori_fn, pOri2, 2, minO, maxO);
// 			if (idx == match_buffer[i][1]){
// 				if (mn == i) mn++;
// 				else{
// 					memcpy(match_buffer[mn], match_buffer[i], 2 * sizeof(int));
// 					mn++;
// 				}
// 			}
// 		}
		LogPrint(0, "Find nearest neighbors in negative direction...");
		int* l_in_r = new int[num_rinl];
		int num_linr = GetBestConditionMatch(l_in_r, pDes1, pLoc1, pScl1, pOri1, num_rinl,
			m_flann_descriptors_index[1], desdistmax, ratiomax,
			er_fn, pLoc2, 2, M, distmax, num_fn, true,
			chek_scale_fn, pScl2, scl_buf_space, minS, maxS,
			chek_ori_fn, pOri2, ori_buf_space, minO, maxO);

		if (num_linr > 1){
			num_linr = 0;
			for (i = 0; i < num_rinl; i++){
				if (match_buffer[i][1] == l_in_r[i]){
					if (num_linr == i) num_linr++;
					else{
						memcpy(match_buffer[num_linr], match_buffer[i], 2 * sizeof(int));
						num_linr++;
					}
				}
			}
		}

		LogPrint(0, "final correspond num = %d", num_linr);
		delete	pDes1;
		delete	pLoc1;
		if (pScl1) delete pScl1;
		delete	l_in_r;
		return num_linr;
	}
// 	int		GetPointMatch( 
// 		BYTE* query_descrpiptors,  float* query_loc, float* query_scale, float* query_ori,
// 		flann_index_t index, float desdistmax, float ratiomax,
// 		transform_err_fn* er_fn, float* loc, int loc_buf_space, void** M, float* distmax, int num_fn, bool bInvert,
// 		check_feature_fn chek_scale_fn, float* scl,int scl_buf_space, float minS, float maxS,
// 		check_feature_fn chek_ori_fn, float* ori,int ori_buf_space, float minO, float maxO
// 		){
// 		enum{ CANDIDATE_POINT_NUM = 20};
// 		int indices[CANDIDATE_POINT_NUM];	float dist[CANDIDATE_POINT_NUM];
// 		int n = flann_radius_search_byte(index, query_descrpiptors, indices, dist, CANDIDATE_POINT_NUM, desdistmax, &m_flann_parameters);
// 		if (n < 1) return -1;
// 		int cnt = 0; int* pIndices = indices;	float* pDist = dist;
// 		int i, j;
// 		for ( i = 0; i < n; i++){
// 			if (er_fn){
// 				float* pLoc = loc + loc_buf_space*indices[i];
// 				float* pLoc1, *pLoc2;
// 				if (bInvert) {
// 					pLoc1 = pLoc;
// 					pLoc2 = query_loc;
// 				}
// 				else{
// 					pLoc2 = pLoc;
// 					pLoc1 = query_loc;
// 				}
// 				for (j = 0; j < num_fn; j++){
// 					if (!er_fn[j](pLoc1, pLoc2, M[j], distmax[j])) break;
// 				}
// 				if (j < num_fn) continue;
// 			}			
// 			if (chek_scale_fn && !chek_scale_fn(*query_scale, *(scl + scl_buf_space*indices[i]), minS, maxS)) continue;
// 			if (chek_ori_fn && !chek_ori_fn(*query_ori, *(ori + ori_buf_space*indices[i]), minO, maxO)) continue;
// 			cnt++;
// 			*pIndices++ = indices[i];
// 			*pDist++ = dist[i];
// 		}
// 		if (cnt<1) return -1;
// 		if (cnt>1){
// 			if (dist[0] > dist[1]*ratiomax) return -1;
// 		}
// 		return indices[0];
// 	}
	int		GetBestConditionMatch(
		int* find_index,
		BYTE* query_descrpiptors, float* query_loc, float* query_scale, float* query_ori, int num,
		void* index, float desdistmax, float ratiomax,
		transform_err_fn* er_fn, float* loc, int loc_buf_space, void** M, float* distmax, int num_fn, bool bInvert,
		check_feature_fn chek_scale_fn, float* scl,int scl_buf_space, float minS, float maxS,
		check_feature_fn chek_ori_fn, float* ori,int ori_buf_space, float minO, float maxO
		)
	{
		int nn = 10;
		int* result = (int*)malloc(num*nn*sizeof(int));
		float* dists = (float*)malloc(num*nn*sizeof(float));
		int i;	float*	pDist = dists;
		for (i = 0; i < num*nn; i++) *pDist++ = -1;
		printf("<1>find nearest neighbors(%d) of feature...", nn);
		FindDesNearestNeighborsIndex(index, query_descrpiptors, num, result, dists, nn);
		printf("done.");

		int cnt_match = 0;
		
		int*	pResult = result;
		pDist = dists;
		float*		pQuery_loc = query_loc;
		float*		pQuery_scale = query_scale;
		float*		pQurey_ori = query_ori;
		int step = 1;
		int step_len = num / 5;	if (step_len <= 0) step_len = 1;
		printf("<2>find best neighbor with guided info");
		for (i = 0; i < num; i++, pDist += nn, pResult += nn, find_index++, pQuery_loc += loc_buf_space, pQuery_scale += scl_buf_space, pQurey_ori += ori_buf_space){
			*find_index = -1;
			if (!((i + 1) % step_len)) {
				printf("...%d", step * 20);	step++;
			}
			int cnt_candi = 0;
			for (int idx = 0; idx < nn; idx++){
				if (cnt_candi >= 2 || pDist[idx] <0 || pResult[idx] < 0 || pDist[idx]>desdistmax) break;
				if (er_fn){
					float* pLoc = loc + loc_buf_space*pResult[idx];
					float* pLoc1, *pLoc2;
					if (bInvert) {
						pLoc2 = pLoc;
						pLoc1 = pQuery_loc;
					}
					else{
						pLoc1 = pLoc;
						pLoc2 = pQuery_loc;						
					}
					int j;
					for ( j = 0; j < num_fn; j++){
						if (!er_fn[j](pLoc1, pLoc2, M[j], distmax[j])) break;
					}
					if (j < num_fn) continue;
				}			
				if (chek_scale_fn){
					if (bInvert){
						if (!chek_scale_fn(*pQuery_scale, *(scl + scl_buf_space*pResult[idx]), minS, maxS)) continue;
					}
					else if (!chek_scale_fn(*(scl + scl_buf_space*pResult[idx]), *pQuery_scale, minS, maxS)) continue;
				}
				if (chek_ori_fn){
					if (bInvert){
						if (!chek_ori_fn(*pQurey_ori, *(ori + ori_buf_space*pResult[idx]), minO, maxO))	continue;
					}
					else if (!chek_ori_fn(*(ori + ori_buf_space*pResult[idx]), *pQurey_ori, minO, maxO))	continue;
				}
				pResult[cnt_candi] = pResult[idx];
				pDist[cnt_candi] = pDist[idx];
				cnt_candi++;
			}
			if (cnt_candi<1) continue;
			if (cnt_candi>1){
				if (pDist[0] > pDist[1] * ratiomax) continue;
			}

			*find_index = *pResult;
			cnt_match++;
		}
		if (step < 6) printf("...100");
		free(result);
		free(dists);
		printf("...done.\n");
		return cnt_match;
	}
	int		GetBestMatch(int* find_index,void* index,BYTE* query_descrpiptors,int num,float distmax,float ratiomax){
		int nn = 2;
		int* result = (int*)malloc(num*nn*sizeof(int));
		float* dists = (float*)malloc(num*nn*sizeof(float));
		FindDesNearestNeighborsIndex(index, query_descrpiptors, num, result, dists, nn);

		int cnt = 0;
		float*	pDist = dists;
		int*	pResult = result;
		for (int i = 0; i < num; i++, pDist += nn, pResult += nn, find_index++){
			*find_index = -1;
			if (*pDist>distmax) continue;
			if (*pDist > *(pDist + 1)*ratiomax) continue;
			*find_index = *pResult;
			cnt++;
		}

		free(result);
		free(dists);
		return cnt;
	}
};

#ifndef _DISTIDX
#define _DISTIDX
typedef struct tagDistIdx{
	int idx;
	float dist;
}DistIdx;
#endif
int		CSiftMatchCPU::GetNeighborSiftMatch(int match_buffer[][2], int maxNeighborNum /* = 2 */, float distmax /* = 0.7f */, int mutual_best_match /* = 1 */){
	Prepare(&distmax, NULL);
	int nn = maxNeighborNum;
	int num = m_sift_num[1];
	int* result = (int*)malloc(num*nn*sizeof(int));
	float* dists = (float*)malloc(num*nn*sizeof(float));
	int i;	float*	pDist = dists;
	for (i = 0; i < num*nn; i++) *pDist++ = -1;

	FindDesNearestNeighborsIndex(m_flann_descriptors_index[0], m_descriptors[1].data(), num, result, dists, nn);

	int cnt = 0;		
	pDist = dists;
	int*	pResult = result;
	if (!mutual_best_match){		
		for ( i = 0; i < num; i++, pDist += nn, pResult += nn){
			for (int j = 0; j < nn; j++){
				if (pDist[j]<0 || pDist[j]>distmax) break;
				(*match_buffer)[1] = i;
				(*match_buffer)[0] = pResult[j];
				cnt++;	match_buffer++;
			}
		}
	}
	else{
		num = m_sift_num[0];
		DistIdx* neighbor_list = (DistIdx*)malloc(num*(nn + 1)*sizeof(DistIdx));
		DistIdx* pNeighbor = neighbor_list;
		float dist_init = distmax + 100;
		for (i = 0; i < num*(nn + 1); i++) {
			pNeighbor->dist = dist_init;
			pNeighbor->idx = -1;
			pNeighbor++;
		}

		for (i = 0; i < m_sift_num[1]; i++, pDist += nn, pResult += nn){		
			for (int j = 0; j < nn; j++){
				if (pDist[j]<0 || pDist[j]>distmax) break;
				pNeighbor = neighbor_list + (pResult[j] + 1)*(nn + 1) - 2;
				float dist_cur = pDist[j];
				int k;
				for (k = 0; k < nn; k++, pNeighbor--){
					if (pNeighbor->dist > dist_cur){
						if (pNeighbor->idx >= 0){
							memcpy(pNeighbor + 1, pNeighbor, sizeof(DistIdx));
						}
					}
					else {
						(pNeighbor + 1)->dist = dist_cur;
						(pNeighbor + 1)->idx = i;
						break;
					}
				}
				if (k == nn){
					(pNeighbor + 1)->dist = dist_cur;
					(pNeighbor + 1)->idx = i;
				}
			}
		}
		pNeighbor = neighbor_list;
		for (i = 0; i < m_sift_num[0]; i++, pNeighbor += nn + 1){
			for (int j = 0; j < nn; j++){
				if ((pNeighbor + j)->idx < 0) break;
				(*match_buffer)[0] = i;
				(*match_buffer)[1] = (pNeighbor + j)->idx;
				match_buffer++;
				cnt++;
			}
		}
		free(neighbor_list);
	}
	

	free(result);
	free(dists);
	LogPrint(0, "Get neighbor[%d] sift match num = %d", maxNeighborNum, cnt);
	return cnt;
}

#ifdef _SIFTGPU
class CSiftMatchGPU : public CSiftMatch
{
public:
	CSiftMatchGPU(){
		m_matcher = g_siftgpu_dll_register.pCreateNewSiftMatchGPU(m_siftmatch_max_num);
	}
	virtual ~CSiftMatchGPU(){
		if (m_matcher) {
			m_matcher->DestroyContextGL();
			delete m_matcher;
		}
	}
	virtual void	SetMaxSift(int max_sift){
		m_siftmatch_max_num = max_sift;
		m_matcher->SetMaxSift(m_siftmatch_max_num);
	}
	virtual void	ParseParam(int argc, char **argv)
	{
		m_matcher->SetDeviceParam(argc, argv);
	}
	virtual bool	_InitEnvi()
	{
		if (!m_matcher) return false;
		if (!GlobalParam::g_siftgpu_initial) {
			if (m_matcher->CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED){
				LogPrint(ERR_ATMCH_FLAG, "Fail to initialize SIFTGPU!");
				return false;
			}
			GlobalParam::g_siftgpu_initial = true;
		}
		else
			m_matcher->VerifyContextGL();
		return true;
	}
	virtual void	SetDescriptors(int index, int num, const unsigned char * descriptors, int id /* = -1 */)
	{
		m_matcher->SetDescriptors(index, num, descriptors, id);
		
		m_sift_num[index] = num;
	}
	virtual void	SetDescriptors(int index, int num, const float* descriptors, int id){
		m_matcher->SetDescriptors(index, num, descriptors, id);
		
		m_sift_num[index] = num;
	}
	virtual void	SetLocation(int index, const float* locations, int gap)
	{
		m_matcher->SetFeautreLocation(index, locations, gap);
	}
	virtual void	SetFeature(int index, const float* feature, int gap)
	{
		m_feature[index].resize(m_sift_num[index] * 2);
		const float* pSrc = feature;
		float* pDst = &m_feature[index][0];
		for (int i = 0; i < m_sift_num[index]; i++){
			memcpy(pDst, pSrc, 2 * sizeof(float));
			pDst += 2;
			pSrc += 2 + gap;
		}
	}
	virtual int		GetSiftMatch(int match_buffer[][2], float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
	{
		int sz =  m_matcher->GetSiftMatch(m_sift_num[0] < m_sift_num[1] ? m_sift_num[0] : m_sift_num[1], match_buffer, distmax, ratiomax, mutual_best_match);
		LogPrint(0, "Get SIFTGPU Correspond Num= %d", sz);
		return sz;
	}
	virtual int		GetNeighborSiftMatch(int match_buffer[][2], int maxNeighborNum /* = 2 */, float distmax /* = 0.7f */, int mutual_best_match /* = 1 */)
	{
		return 0;
	}
	virtual int		GetGuidedSiftMatch(int match_buffer[][2], float H[3][3], float F[3][3], float distmax, float ratiomax, float hdistmax, float fdistmax, int mutual_best_match)
	{
		int sz = m_matcher->GetGuidedSiftMatch(m_sift_num[0] < m_sift_num[1] ? m_sift_num[0] : m_sift_num[1], match_buffer, H, F, distmax, ratiomax, hdistmax, fdistmax, mutual_best_match);
		LogPrint(0, "Get guided SIFTGPU Correspond Num= %d", sz);
		return sz;
	}
	virtual int GetXformSiftMatch(
		int match_buffer[][2],
		transform_err_fn er_fn, void* xform, float errdistmax,
		float distmax,
		float ratiomax,
		int mutual_best_match
		)
	{
		return 0;
	}
	virtual int		GetLimitedSiftMatch(int match_buffer[][2], float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
	{
		return 0;
	}
	virtual int		GetLimitedGuidedSiftMatch(int match_buffer[][2], float H[3][3], float F[3][3], float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax, float ratiomax, float hdistmax, float fdistmax, int mutual_best_match)
	{
		return 0;
	}
	virtual int		GetLimitedXformSiftMatch(int match_buffer[][2], transform_err_fn er_fn, void* xform, float errdistmax, float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax, float ratiomax, int mutual_best_match)
	{
		return 0;
	}
private:
	SiftMatchGPU*			m_matcher;
	vector<float>			m_feature[2];
	int						m_sift_num[2];
};
#endif

bool		CheckLimitedScale(float ratiomin_scale, float ratiomax_scale)
{
	return (ratiomin_scale >= 0 && ratiomax_scale > 0);
}
bool		CheckLimitedOrientation(float diffmin_ori, float diffmax_ori){
	return (diffmax_ori <= PI && diffmin_ori >= 0);
}
#define MATCH_MODEL_NOTHING		0
#define MATCH_MODEL_DESCRIP		1
#define MATCH_MODEL_LOCATION	2
#define MATCH_MODEL_FEATURE		4

CSiftMatch::CSiftMatch(){
	m_siftmatch_max_num = 8192;//
	m_sift_matcher = NULL;
	m_match_condition[0] = m_match_condition[1] = MATCH_MODEL_NOTHING;
}

CSiftMatch::~CSiftMatch(){
	if (m_sift_matcher) delete m_sift_matcher;
}

void* CSiftMatch::operator new (size_t  size){
	void * p = malloc(size);
	if (p == 0)
	{
		const std::bad_alloc ba;
		throw ba;
	}
	return p;
}

bool CSiftMatch::InitEnvi(bool bGPU){
	if (m_sift_matcher)  {
		delete m_sift_matcher; m_sift_matcher = NULL;
	}
	m_match_condition[0] = m_match_condition[1] = MATCH_MODEL_NOTHING;
#ifdef _SIFTGPU
	if (bGPU) {
		m_sift_matcher = new CSiftMatchGPU;
		if (m_sift_matcher){
			if (m_sift_matcher->_InitEnvi()){
				LogPrint(0, "[USE SIFTMATCHGPU]");
				return true;
			}
			delete m_sift_matcher; m_sift_matcher = NULL;
		}
		LogPrint(0, "SIFTGPU is not supported.Check if exist SIFTGPU library.");
	}
#endif
	LogPrint(0, "[USE SIFTMATCHCPU]");
	m_sift_matcher = new CSiftMatchCPU;
	return m_sift_matcher ? m_sift_matcher->_InitEnvi() : false;
}

void	CSiftMatch::SetMaxSift(int max_sift)
{
	if (!m_sift_matcher&&!InitEnvi(GlobalParam::g_sift_match_gpu)) return;
	m_sift_matcher->SetMaxSift(max_sift);
}

void	CSiftMatch::ParseParam(int argc, char **argv)
{
	if (!m_sift_matcher&&!InitEnvi(GlobalParam::g_sift_match_gpu)) return;
	m_sift_matcher->ParseParam(argc,argv);
}

void	CSiftMatch::SetDescriptors(int index, int num, const unsigned char * descriptors, int id /* = -1 */)
{
	if (!m_sift_matcher&&!InitEnvi(GlobalParam::g_sift_match_gpu)) return;
	if (index > 1) index = 1; else
	if (index < 0) index = 0;
	m_match_condition[index] = MATCH_MODEL_DESCRIP;
	m_sift_matcher->SetDescriptors(index, num, descriptors, id);
}

void	CSiftMatch::SetDescriptors(int index, int num, const float* descriptors, int id /* = -1 */)
{
	if (!m_sift_matcher&&!InitEnvi(GlobalParam::g_sift_match_gpu)) return;
	if (index > 1) index = 1; else
	if (index < 0) index = 0;
	m_match_condition[index] = MATCH_MODEL_DESCRIP;
	m_sift_matcher->SetDescriptors(index, num, descriptors, id);
}

void	CSiftMatch::SetLocation(int index, const float* locations, int gap /* = 0 */)
{
	if (!m_sift_matcher&&!InitEnvi(GlobalParam::g_sift_match_gpu)) return;
	if (index > 1) index = 1; else
	if (index < 0) index = 0;
	m_match_condition[index] |= MATCH_MODEL_LOCATION;
	m_sift_matcher->SetLocation(index, locations, gap);
}

void	CSiftMatch::SetFeature(int index, const float* feature, int gap /* = 0 */)
{
	if (!m_sift_matcher&&!InitEnvi(GlobalParam::g_sift_match_gpu)) return;
	if (index > 1) index = 1; else
	if (index < 0) index = 0;
	m_match_condition[index] |= MATCH_MODEL_FEATURE;
	m_sift_matcher->SetFeature(index, feature, gap);
}

int		CSiftMatch::GetSiftMatch(int match_buffer[][2], float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	return m_sift_matcher->GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
}

int		CSiftMatch::GetNeighborSiftMatch(int match_buffer[][2], int maxNeighborNum /* = 2 */, float distmax /* = 0.7f */, int mutual_best_match /* = 1 */)
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	return m_sift_matcher->GetNeighborSiftMatch(match_buffer, maxNeighborNum, distmax, mutual_best_match);
}

int		CSiftMatch::GetGuidedSiftMatch(int match_buffer[][2], float H[3][3], float F[3][3], float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, float hdistmax /* = 32 */, float fdistmax /* = 16 */, int mutual_best_match /* = 1 */)
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	if (H == NULL && F == NULL || (!(m_match_condition[0] & MATCH_MODEL_LOCATION) || !(m_match_condition[1] & MATCH_MODEL_LOCATION)))
	{
		return m_sift_matcher->GetSiftMatch( match_buffer, distmax, ratiomax, mutual_best_match);
	}
	else
	{
		float Z[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }, ti = (1.0e+20F);

		return m_sift_matcher->GetGuidedSiftMatch( match_buffer, H ? H : Z, F ? F : Z,
			distmax, ratiomax, H ? hdistmax : ti, F ? fdistmax : ti, mutual_best_match);
	}
}

int		CSiftMatch::GetXformSiftMatch( int match_buffer[][2], transform_err_fn er_fn, void* xform, float errdistmax, float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */ )
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	if (!er_fn || (!(m_match_condition[0] & MATCH_MODEL_LOCATION) || !(m_match_condition[1] & MATCH_MODEL_LOCATION)))
	{
		return m_sift_matcher->GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
	}
	else
	{
		
		return m_sift_matcher->GetXformSiftMatch(match_buffer, er_fn, xform, errdistmax, distmax, ratiomax , mutual_best_match);
	}
}

int		CSiftMatch::GetLimitedSiftMatch(int match_buffer[][2], float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	bool bScale = CheckLimitedScale(ratiomin_scale, ratiomax_scale);
	bool bOri = CheckLimitedOrientation(diffmin_ori, diffmax_ori);
	if ( !bScale && !bOri || (!(m_match_condition[0] & MATCH_MODEL_FEATURE) || !(m_match_condition[1] & MATCH_MODEL_FEATURE))){
		return m_sift_matcher->GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
	}
	else
	{
		float rmin_scale = 1.0e-6f,rmax_scale = 1.0e+20f;
		float rmin_ori = -0.01f, rmax_ori = PI + 0.01f;
		diffmin_ori -= 0.01f; diffmax_ori += 0.01f;
		return m_sift_matcher->GetLimitedSiftMatch(match_buffer, bScale ? ratiomin_scale : rmin_scale, bScale ? ratiomax_scale : rmax_scale, bOri ? diffmin_ori : rmin_ori, bOri ? diffmax_ori : rmax_ori , distmax, ratiomax, mutual_best_match);
	}
	
}

int		CSiftMatch::GetLimitedGuidedSiftMatch( int match_buffer[][2], float H[3][3], float F[3][3], float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, float hdistmax /* = 32 */, float fdistmax /* = 16 */, int mutual_best_match /* = 1 */ )
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	int model_loc_fea = MATCH_MODEL_DESCRIP | MATCH_MODEL_LOCATION | MATCH_MODEL_FEATURE;
	bool bLoc = H || F;
	bool bScale = CheckLimitedScale(ratiomin_scale, ratiomax_scale);
	bool bOri = CheckLimitedOrientation(diffmin_ori, diffmax_ori);
	float rmin_scale = 1.0e-6f, rmax_scale = 1.0e+20f;
	float rmin_ori = -0.01f, rmax_ori = PI + 0.01f;
	diffmin_ori -= 0.01f; diffmax_ori += 0.01f;
	float Z[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	if ((m_match_condition[0] == model_loc_fea) && (m_match_condition[1] == model_loc_fea) && bLoc && (bScale || bOri)){		
		return m_sift_matcher->GetLimitedGuidedSiftMatch(match_buffer, H ? H : Z, F ? F : Z,
			bScale ? ratiomin_scale : rmin_scale, bScale ? ratiomax_scale : rmax_scale, bOri ? diffmin_ori : rmin_ori, bOri ? diffmax_ori : rmax_ori,
			distmax, ratiomax, H ? hdistmax : rmax_scale, F ? fdistmax : rmax_scale, mutual_best_match);
	}
	else if ((m_match_condition[0] & MATCH_MODEL_FEATURE) && (m_match_condition[1] & MATCH_MODEL_FEATURE) && (bScale || bOri) ){
		return m_sift_matcher->GetLimitedSiftMatch(match_buffer, bScale ? ratiomin_scale : rmin_scale, bScale ? ratiomax_scale : rmax_scale, bOri ? diffmin_ori : rmin_ori, bOri ? diffmax_ori : rmax_ori, distmax, ratiomax, mutual_best_match);
	}
	else if ((m_match_condition[0] & MATCH_MODEL_LOCATION) && (m_match_condition[1] & MATCH_MODEL_LOCATION)&& bLoc ){
		return m_sift_matcher->GetGuidedSiftMatch(match_buffer, H ? H : Z, F ? F : Z,
			distmax, ratiomax, H ? hdistmax : rmax_scale, F ? fdistmax : rmax_scale, mutual_best_match);
	}
	else{
		return m_sift_matcher->GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
	}
}

int		CSiftMatch::GetLimitedXformSiftMatch(int match_buffer[][2], transform_err_fn er_fn, void* xform, float errdistmax, float ratiomin_scale, float ratiomax_scale, float diffmin_ori, float diffmax_ori, float distmax /* = 0.7 */, float ratiomax /* = 0.8 */, int mutual_best_match /* = 1 */)
{
	if (!m_sift_matcher) return 0;
	if (!(m_match_condition[0] & MATCH_MODEL_DESCRIP) || !(m_match_condition[1] & MATCH_MODEL_DESCRIP)) return 0;
	int model_loc_fea = MATCH_MODEL_DESCRIP | MATCH_MODEL_LOCATION | MATCH_MODEL_FEATURE;
	bool bScale = CheckLimitedScale(ratiomin_scale, ratiomax_scale);
	bool bOri = CheckLimitedOrientation(diffmin_ori, diffmax_ori);
	float rmin_scale = 1.0e-6f, rmax_scale = 1.0e+20f;
	float rmin_ori = -0.01f, rmax_ori = PI + 0.01f;
	diffmin_ori -= 0.01f; diffmax_ori += 0.01f;
	if ((m_match_condition[0] == model_loc_fea) && (m_match_condition[1] == model_loc_fea) && er_fn && (bScale || bOri)){
		return m_sift_matcher->GetLimitedXformSiftMatch(match_buffer, er_fn,xform,errdistmax,
			bScale ? ratiomin_scale : rmin_scale, bScale ? ratiomax_scale : rmax_scale, bOri ? diffmin_ori : rmin_ori, bOri ? diffmax_ori : rmax_ori,
			distmax, ratiomax, mutual_best_match);
	}
	else if ((m_match_condition[0] & MATCH_MODEL_FEATURE) && (m_match_condition[1] & MATCH_MODEL_FEATURE) && (bScale || bOri)){
		return m_sift_matcher->GetLimitedSiftMatch(match_buffer, bScale ? ratiomin_scale : rmin_scale, bScale ? ratiomax_scale : rmax_scale, bOri ? diffmin_ori : rmin_ori, bOri ? diffmax_ori : rmax_ori, distmax, ratiomax, mutual_best_match);
	}
	else if ((m_match_condition[0] & MATCH_MODEL_LOCATION) && (m_match_condition[1] & MATCH_MODEL_LOCATION) && er_fn){
		return m_sift_matcher->GetXformSiftMatch(match_buffer, er_fn, xform, errdistmax, distmax, ratiomax, mutual_best_match);
	}
	else{
		return m_sift_matcher->GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
	}
}

int		GetComboSiftMatch(
	const float* locationsL, const BYTE* desL, int ptNumL,
	const float* locationsR, const BYTE* desR, int ptNumR,
	int match_buffer[][2],
	bool bUseH , bool bUseF,
	float distmax ,
	float ratiomax,
	int mutual_best_match,
	bool bGPU
	)
{
	CSiftMatch matcher;
	if (!matcher.InitEnvi(bGPU)) return false;

	matcher.SetDescriptors(0, ptNumL, desL);
	matcher.SetDescriptors(1, ptNumR, desR);

	if (!locationsL || !locationsR || (!bUseF&&!bUseH)){
		return matcher.GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
	}

	LogPrint(0, "[step 1]GetSiftMatch...");
	int sift_num = matcher.GetSiftMatch(match_buffer, distmax, ratiomax, 1);
	LogPrint(0, "Get sift correspond points num = %d", sift_num);

	if ( sift_num < 4 || (bUseF&&sift_num<8) ){
		LogPrint(0, "too less match points.so terminate combo match.");
		return sift_num;
	}
	LogPrint(0, "[step 2]Calculate %s %s ...", bUseH ? "H" : "", bUseF ? "F" : "");
	float H[3][3], F[3][3];
	int nBufferSpace[4] = { 2, 2, 2, 2 };
	
	char * argv[] = { "method", "ransacA", "matrixCalc", "inliers", "ProbBadSupp", "0.01", "DistanceThreshold","16" };
	int argc = sizeof(argv) / sizeof(char*);
	int h_left_num = 0;
	if(bUseH) h_left_num = GeoTransform::EstimateXformMatrix(locationsL, locationsL + 1, locationsR, locationsR + 1, nBufferSpace, match_buffer, sift_num, NULL, H[0],
#ifdef GEOTRANSFORM_STDARG
		1, "DistanceThreshold", 16.0);
#else
		argc, argv);
#endif
// 	if (h_left_num){
// 		LogPrint(0, "h[3][3] = [%5.2lf %5.2lf %5.2lf; %5.2lf %5.2lf %5.2lf; %5.2lf %5.2lf %5.2lf]",
// 			H[0][0], H[0][1], H[0][2], H[1][0], H[1][1], H[1][2], H[2][0], H[2][1], H[2][2]);
// 	}
	argv[argc - 1] = "8";
	int f_left_num = 0;
	if(bUseF) f_left_num = GeoTransform::EstimateFundamentalMatrix(locationsL, locationsL + 1, locationsR, locationsR + 1, nBufferSpace, match_buffer, sift_num, NULL, F[0], 
#ifdef GEOTRANSFORM_STDARG
		1, "DistanceThreshold", 8.0);
#else
		argc, argv);
#endif
// 	if (f_left_num){
// 		LogPrint(0, "f[3][3] = [%5.2lf %5.2lf %5.2lf; %5.2lf %5.2lf %5.2lf; %5.2lf %5.2lf %5.2lf]",
// 			F[0][0], F[0][1], F[0][2], F[1][0], F[1][1], F[1][2], F[2][0], F[2][1], F[2][2]);
// 	}

	if (!h_left_num&&!f_left_num){
		LogPrint(0, "fail to get h or f.so terminate combo match.");
		return sift_num;
	}

	LogPrint(0, "[step 3]GetGuidedSiftMatch...");

	matcher.SetLocation(0, locationsL);
	matcher.SetLocation(1, locationsR);

	return matcher.GetGuidedSiftMatch(match_buffer, h_left_num ? H : NULL, f_left_num ? F : NULL, distmax, ratiomax, 32.0f, 16.0f, mutual_best_match);
}

int		GetSimilaritySiftMatch( 
	const float* locationsL, const BYTE* desL, const float* featureL, int ptNumL,
	const float* locationsR, const BYTE* desR, const float* featureR, int ptNumR,
	int match_buffer[][2],
	float distmax ,
	float ratiomax ,
	int mutual_best_match,
	bool bGPU
	)
{
	CSiftMatch matcher;
	if (!matcher.InitEnvi(bGPU)) return false;

	matcher.SetDescriptors(0, ptNumL, desL);
	matcher.SetDescriptors(1, ptNumR, desR);

	if (!locationsL || !locationsR){
		return matcher.GetSiftMatch(match_buffer, distmax, ratiomax, mutual_best_match);
	}

	LogPrint(0, "[step 1]GetSiftMatch...");
	int sift_num = matcher.GetSiftMatch(match_buffer, distmax, ratiomax, 1);
	LogPrint(0, "Get sift correspond points num = %d", sift_num);

	if (sift_num < 8){
		LogPrint(0, "too less match points.so terminate combo match.");
		return sift_num;
	}
	LogPrint(0, "[step 2]Calculate Similarity transform...");
	float H[3][3];
	int nBufferSpace[4] = { 2, 2, 2, 2 };
	char * argv[] = { "xform", "nonreflective similarity", "DistanceThreshold", "16","ProbBadSupp","0.6" };
	int argc = sizeof(argv) / sizeof(char*);
	int h_left_num = GeoTransform::EstimateXformMatrix(locationsL, locationsL + 1, locationsR, locationsR + 1, nBufferSpace, match_buffer, sift_num, NULL, H[0],
#ifdef GEOTRANSFORM_STDARG
		2,"xform","nonreflective similarity", "DistanceThreshold", 16.0);
#else
		argc, argv);
#endif
	if (!h_left_num){
		LogPrint(0, "fail to get h.so terminate combo match.");
		return sift_num;
	}
	
	LogPrint(0, "[step 3]GetGuidedSiftMatch...");
	float scale = (float)sqrt((H[0][0] * H[0][0] + H[0][1] * H[0][1] + H[1][0] * H[1][0] + H[1][1] * H[1][1]) / 2);
	float kappa = 0;
	if (fabs(H[1][0]) < 1e-6) kappa = 0;
	else if (fabs(H[0][0]) < 1e-6) kappa = PI/2;
	else {
		kappa = (float)fabs(atan2(H[1][0] , H[0][0]));
	}
	float ori_min = kappa - PI / 6;	if (ori_min < 0) ori_min = 0;
	float ori_max = kappa + PI / 6;	if (ori_max > PI) ori_min = PI;
	matcher.SetLocation(0, locationsL);
	matcher.SetLocation(1, locationsR);
	matcher.SetFeature(0, featureL);
	matcher.SetFeature(1, featureR);

	return matcher.GetLimitedGuidedSiftMatch(match_buffer, H, NULL, scale / 1.414f, scale*1.414f, ori_min, ori_max, distmax, ratiomax, 32.0f, 0, mutual_best_match);
}

void	FreeSiftPtsMem(void* pt)
{
	if (pt) delete[] pt;
}

#include "algorithms/Pretreatment.h"

int		ExtractSiftFromBuffer(const BYTE* pImg, int nCols, int nRows, float** loc0, float** fea0, BYTE** des0,bool bGPU){
	int nColNum, nRowNum;

	CSift sift;
	if ( !sift.InitEnvi(bGPU) ) return 0;
	int* rc_split = Pretreatment::split_image(nColNum, nRowNum, 0, 0, nCols, nRows, 8, (int)(GlobalParam::g_extract_buffer_size), Pretreatment::RECT_FIRST);

	BYTE* pTmpBuf = new BYTE[GlobalParam::g_extract_buffer_size];

	int *pSplitRc = rc_split;
	
	int ptSum = 0;
	int sum = nColNum*nRowNum;
	float* loc = NULL, *fea = NULL; BYTE* des = NULL;
	for (int i = 0; i <sum ; i++, pSplitRc += 4){
		NORMALISE_IMAGE_SIDE(pSplitRc[2]);		NORMALISE_IMAGE_SIDE(pSplitRc[3]);
		Pretreatment::ReadBlock4Buf(pImg, nCols, nRows, pSplitRc[0], pSplitRc[1], pTmpBuf, pSplitRc[2], pSplitRc[3]);
		printf( "extract sift block[%d/%d]...\n", i + 1, sum);
		if (sift.RunSIFT(pTmpBuf, pSplitRc[2], pSplitRc[3])) {
			int num = sift.GetFeatureNum();
			if (num < 1) continue;
			int last_num = ptSum;
			ptSum += num;
			float* loc_last = new float[ptSum * 2];	memcpy(loc_last, loc, sizeof(float)* 2 * last_num);
			if(loc) delete loc;	loc = loc_last;
			float* fea_last = new float[ptSum * 2];	memcpy(fea_last, fea, sizeof(float)* 2 * last_num);
			if(fea) delete fea;	fea = fea_last;
			BYTE*  des_last = new BYTE[ptSum * 128];	memcpy(des_last, des, sizeof(BYTE)* 128 * last_num);
			if(des) delete des;	des = des_last;
						
			sift.GetFeatureVector(loc + last_num * 2, fea + 2 * last_num, des + 128 * last_num);
			float* pLoc = loc + last_num * 2;
			for (int j = 0; j < num; j++){
				*pLoc += pSplitRc[0];	pLoc++;
				*pLoc += pSplitRc[1];	pLoc++;
			}
		}
	}
	delete pTmpBuf;
	delete rc_split;
	if (loc0) *loc0 = loc; else delete loc;
	if (fea0) *fea0 = fea; else delete fea;
	if (des0) *des0 = des; else delete des;
	LogPrint(0, "extract sift num = %d", ptSum);
	return ptSum;
}

void		FreeSiftFeature(float* loc, float* fea, BYTE* des){
	if (loc) delete loc;
	if (fea) delete fea;
	if (des) delete des;
}

#ifdef _DEBUG
bool SaveSift(const char* lpstrPathName,float* loc,int sz,int nRows){
	FILE* fp = fopen(lpstrPathName, "w");	if (!fp) return false;
	fprintf(fp, "%d\n",sz);
	for (int i = 0; i < sz; i++ ){
		fprintf(fp, "%5d\t%f\t%f\t%f\n", i, *loc, nRows - 1 - *(loc + 1), *(loc + 1));
		loc += 2;
	}
	fclose(fp);
	return true;
}
#endif

int ascend_cmp_FPT4D(const void* a, const void* b){
	FPT4D* pA = (FPT4D*)a;	FPT4D* pB = (FPT4D*)b;
	float dx = pA->xl - pB->xl;
	if (dx>1e-5){
		return 1;
	}else if (dx<-1e-5){
		return -1;
	}
	else{
		float dy = pA->yl - pB->yl;
		if (dy >= 0) return 1;
		else return -1;
	}
}

FPT4D* GetSiftMatch(int* ptSum, const  BYTE *pImgL, int colsL, int rowsL, const  BYTE *pImgR, int colsR, int rowsR, bool bUseH, bool bUseF, float distmax, float ratiomax, int mutual_best_match, bool bGPU)
{
	*ptSum = 0;
	FPT4D* mpt = NULL;
	CSiftMatch matcher;

	float *locL,*locR; 
	float *feaL,*feaR;
	BYTE *desL,*desR;
	LogPrint(0, "[step 1]Extract sift feature from left image...");	
	int ptNumL = ExtractSiftFromBuffer(pImgL, colsL, rowsL, &locL, &feaL, &desL,bGPU);
	delete feaL;
// #ifdef _DEBUG
// 	SaveSift("D:/Data/shanxi/ws/pyrml.txt", (float*)locL.data(), ptNumL,rowsL);
// #endif
	if (ptNumL < 3) goto loop1;
	LogPrint(0, "[step 2]Extract sift feature from right image...");
	int ptNumR = ExtractSiftFromBuffer(pImgR, colsR, rowsR, &locR, &feaR, &desR,bGPU);
	delete feaR;
// #ifdef _DEBUG
// 	SaveSift("D:/Data/shanxi/ws/pyrmr.txt", (float*)locR.data(), ptNumR,rowsR);
// #endif
	if (ptNumR < 3) goto loop2;

	LogPrint(0, "[step 3]GetSiftMatch...");
	
	if (!matcher.InitEnvi(bGPU)) goto loop2;
	
	int(*match_buf)[2] = new int[ptNumL>ptNumR ? ptNumL : ptNumR][2];
	matcher.SetDescriptors(0, ptNumL, desL);
	matcher.SetDescriptors(1, ptNumR, desR);
		
	int match_num = matcher.GetSiftMatch(match_buf, distmax, ratiomax, mutual_best_match);

	if (match_num>0 && (bUseF || bUseH) ){
		LogPrint(0, "[step 4]Calculate %s %s ...", bUseH ? "H" : "", bUseF ? "F": "");
		float H[3][3], F[3][3];
		int nBufferSpace[4] = { 2, 2, 2, 2 };
		int h_left_num = 0;	int f_left_num = 0;
		
		char * argv[] = { "method", "ransacA","DistanceThreshold", "16" };
		int argc = sizeof(argv) / sizeof(char*);
		if (bUseH) h_left_num = GeoTransform::EstimateXformMatrix(locL, locL + 1, locR, locR + 1, nBufferSpace, match_buf, match_num, NULL, H[0], 
#ifdef GEOTRANSFORM_STDARG
			1, "DistanceThreshold", 16.0);
#else
			argc, argv);
#endif
		argv[3] = "8";
		
		if (bUseF) f_left_num = GeoTransform::EstimateFundamentalMatrix(locL, locL + 1, locR, locR + 1, nBufferSpace, match_buf, match_num, NULL, F[0],
#ifdef GEOTRANSFORM_STDARG
			1, "DistanceThreshold", 8.0);
#else
			argc, argv);
#endif
		if (h_left_num||f_left_num){
			LogPrint(0, "[step 5]GetGuidedSiftMatch...");

			matcher.SetLocation(0, locL);
			matcher.SetLocation(1, locR);

			match_num = matcher.GetGuidedSiftMatch(match_buf, h_left_num ? H : NULL, f_left_num ? F : NULL, distmax, ratiomax, 32.0f, 16.0f, mutual_best_match);
		}else
			LogPrint(0, "fail to get h or f.");
	}

	if (match_num>0){
		mpt = new FPT4D[match_num];
		FPT4D* pMPT = mpt;
		int i;
		for ( i = 0; i < match_num; i++){
			float* xyl = locL + match_buf[i][0] * 2;
			float* xyr = locR + match_buf[i][1] * 2;
			pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
			pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
			pMPT++;
		}
		qsort(mpt, match_num, sizeof(FPT4D), ascend_cmp_FPT4D);
		float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
		int sz = 1;	pMPT = mpt+1;	FPT4D* pMPT0 = mpt;
		for (i = 1; i < match_num; i++, pMPT++){
			if (fabs(pMPT->xl - pMPT0->xl) < 1e-5 && fabs(pMPT->yl - pMPT0->yl) < 1e-5
				&& fabs(pMPT->xr - pMPT0->xr) < 1e-5 && fabs(pMPT->yr - pMPT0->yr) < 1e-5) continue;
			pMPT0++;	if (pMPT0 != pMPT) memcpy(pMPT0, pMPT, sizeof(FPT4D));
			sz++;
		}
		
		*ptSum = sz;
		LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
	}

	delete[] match_buf;
loop2:
	delete locR;	delete desR;
loop1:
	delete locL;	delete desL;
	return mpt;
}

FPT4D*	 GetNeighborSiftMatch(
	int* ptSum,
	const BYTE *pImgL, int colsL, int rowsL,
	const BYTE *pImgR, int colsR, int rowsR,
	int maxNeighborNum ,
	float distmax,
	int mutual_best_match 
	)
{
	*ptSum = 0;
	FPT4D* mpt = NULL;
	CSiftMatch matcher;
	float *locL, *locR;
	float *feaL, *feaR;
	BYTE *desL, *desR;
	bool bGPU = false;
	LogPrint(0, "[step 1]Extract sift feature from left image...");
	int ptNumL = ExtractSiftFromBuffer(pImgL, colsL, rowsL, &locL, &feaL, &desL, bGPU);
	delete feaL; 
	if (ptNumL < 3) goto loop1;
	LogPrint(0, "[step 2]Extract sift feature from right image...");
	int ptNumR = ExtractSiftFromBuffer(pImgR, colsR, rowsR, &locR, &feaR, &desR, bGPU);
	delete feaR;
	if (ptNumR < 3) goto loop2;

	LogPrint(0, "[step 3]GetNeighborSiftMatch...");
	if (!matcher.InitEnvi(bGPU)) goto loop2;

	int(*match_buf)[2] = new int[(ptNumL>ptNumR ? ptNumL : ptNumR)*maxNeighborNum][2];
	matcher.SetDescriptors(0, ptNumL, desL);
	matcher.SetDescriptors(1, ptNumR, desR);
	int match_num = matcher.GetNeighborSiftMatch(match_buf,maxNeighborNum, distmax, mutual_best_match);

	if (match_num > 0){
		mpt = new FPT4D[match_num];
		FPT4D* pMPT = mpt;
		int i;
		for (i = 0; i < match_num; i++){
			float* xyl = locL + match_buf[i][0] * 2;
			float* xyr = locR + match_buf[i][1] * 2;
			pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
			pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
			pMPT++;
		}
		qsort(mpt, match_num, sizeof(FPT4D), ascend_cmp_FPT4D);
		float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
		int sz = 1;	pMPT = mpt + 1;	FPT4D* pMPT0 = mpt;
		for (i = 1; i < match_num; i++, pMPT++){
			if (fabs(pMPT->xl - pMPT0->xl) < 1e-5 && fabs(pMPT->yl - pMPT0->yl) < 1e-5
				&& fabs(pMPT->xr - pMPT0->xr) < 1e-5 && fabs(pMPT->yr - pMPT0->yr) < 1e-5) continue;
			pMPT0++;	if (pMPT0 != pMPT) memcpy(pMPT0, pMPT, sizeof(FPT4D));
			sz++;
		}

		*ptSum = sz;
		LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
	}

	delete[] match_buf;
loop2:
	delete locR;	delete desR;
loop1:
	delete locL;	delete desL;
	return mpt;
}

FPT4D*	 GetXformSiftMatch(
	int* ptSum,
	const BYTE *pImgL, int colsL, int rowsL,
	const BYTE *pImgR, int colsR, int rowsR,
	transform_err_fn er_fn, void* xform, float errdistmax,
	float distmax ,
	float ratiomax ,
	int mutual_best_match
	)
{
	*ptSum = 0;
	FPT4D* mpt = NULL;
	CSiftMatch matcher;
	float *locL, *locR;
	float *feaL, *feaR;
	BYTE *desL, *desR;
	bool bGPU = false;
	LogPrint(0, "[step 1]Extract sift feature from left image...");
	int ptNumL = ExtractSiftFromBuffer(pImgL, colsL, rowsL, &locL, &feaL, &desL, bGPU);
	delete feaL;
	if (ptNumL < 3) goto loop1;
	LogPrint(0, "[step 2]Extract sift feature from right image...");
	int ptNumR = ExtractSiftFromBuffer(pImgR, colsR, rowsR, &locR, &feaR, &desR, bGPU);
	delete feaR;
	if (ptNumR < 3) goto loop2;

	LogPrint(0, "[step 3]GetNeighborSiftMatch...");
	if (!matcher.InitEnvi(bGPU)) goto loop2;

	int(*match_buf)[2] = new int[ptNumL>ptNumR ? ptNumL : ptNumR][2];
	matcher.SetDescriptors(0, ptNumL, desL);	matcher.SetLocation(0, locL);
	matcher.SetDescriptors(1, ptNumR, desR);	matcher.SetLocation(1, locR);
	int match_num = matcher.GetXformSiftMatch(match_buf, er_fn, xform, errdistmax, distmax, ratiomax, mutual_best_match);

	if (match_num > 0){
		mpt = new FPT4D[match_num];
		FPT4D* pMPT = mpt;
		int i;
		for (i = 0; i < match_num; i++){
			float* xyl = locL + match_buf[i][0] * 2;
			float* xyr = locR + match_buf[i][1] * 2;
			pMPT->xl = *xyl;	pMPT->yl = *(xyl + 1);
			pMPT->xr = *xyr;	pMPT->yr = *(xyr + 1);
			pMPT++;
		}
		qsort(mpt, match_num, sizeof(FPT4D), ascend_cmp_FPT4D);
		float lx0 = -99, ly0 = -99, rx0 = -99, ry0 = -99;
		int sz = 1;	pMPT = mpt + 1;	FPT4D* pMPT0 = mpt;
		for (i = 1; i < match_num; i++, pMPT++){
			if (fabs(pMPT->xl - pMPT0->xl) < 1e-5 && fabs(pMPT->yl - pMPT0->yl) < 1e-5
				&& fabs(pMPT->xr - pMPT0->xr) < 1e-5 && fabs(pMPT->yr - pMPT0->yr) < 1e-5) continue;
			pMPT0++;	if (pMPT0 != pMPT) memcpy(pMPT0, pMPT, sizeof(FPT4D));
			sz++;
		}

		*ptSum = sz;
		LogPrint(0, "Remove repeat correspond points.Final Num= %d", sz);
	}

	delete[] match_buf;
loop2:
	delete locR;	delete desR;
loop1:
	delete locL;	delete desL;
	return mpt;
}