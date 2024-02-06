#pragma once

#include "CSDataDefExport.h"
#include "IFamilyInstance.h"
#include "qCadBaseCmd.h"

#include ".\FindLoop\ISearchLoop.h"


//#define DY_TASK1


const double ARADK = 0.0174532925;
using namespace QGEO;
using namespace gea;

template<class T1, class T2>
struct dmap
{
	std::map<T1, T2>m_mapObj;

	typedef typename std::map<T1, T2> dmap_pair;
	typedef typename dmap_pair::value_type value_type;
	typedef typename dmap_pair::const_iterator const_iterator;
	typedef typename dmap_pair::iterator iterator;

public:
	~dmap()
	{
		clear();
	};

	void insert(T1 id, T2 da, bool ifRepaceExist = false)
	{
		if (ifRepaceExist)
		{
			iterator iter = m_mapObj.find(id);
			if (iter != m_mapObj.end())
			{
				m_mapObj[id] = da;
			}
			else
			{
				m_mapObj.insert(dmap_pair::value_type(id, da));
			}
		}
		else
		{
			m_mapObj.insert(dmap_pair::value_type(id, da));
		}
	};

	bool has(T1 id, T2& da)
	{
		iterator iter = m_mapObj.find(id);
		if (iter != m_mapObj.end())
		{
			da = m_mapObj[id];
			return true;
		}

		return false;
	};

	BOOL getAt(const T1 id, T2& da)
	{
		iterator iter = m_mapObj.find(id);
		if (iter != m_mapObj.end())
		{
			da = m_mapObj[id];
			return TRUE;
		}
		return FALSE;
	};

	void setAt(const T1 id, T2& da)
	{
		iterator iter = m_mapObj.find(id);
		if (iter != m_mapObj.end())
		{
			m_mapObj[id] = da;
		}
		else
		{
			m_mapObj.insert(dmap_pair::value_type(id, da));
		}
	};

	void clear()
	{
		m_mapObj.clear();
	};
};


struct CurveDataAlum
{
	QGEO::geoPoint3d ps;
	QGEO::geoPoint3d pe;
	QGEO::geoPoint3d pm;  // 如果是弧线，需要中间点

	bool isArc;
	DWORD_PTR userData;

	CurveDataAlum()
	{
		ps.Set(0, 0, 0);
		pe.Set(0, 0, 0);
		pm.Set(0, 0, 0);

		isArc = false;
		userData = NULL;
	};

	CurveDataAlum(double x1, double y1, double r, double x2, double y2, DWORD_PTR usd = NULL)
	{
		if (r > 1)
		{
			isArc = true;

			ps.x = x1 + cos(x2 * ARADK) * r;
			ps.y = y1 + sin(x2 * ARADK) * r;

			pe.x = x1 + cos(y2 * ARADK) * r;
			pe.y = y1 + sin(y2 * ARADK) * r;

			double am = x2 + (y2 - x2) * .5;
			pm.x = x1 + cos(am * ARADK) * r;
			pm.y = y1 + sin(am * ARADK) * r;
		}
		else
		{
			isArc = false;

			ps.Set(x1, y1, 0);
			pe.Set(x2, y2, 0);
			pm.Set((x1 + x2) * .5, (y1 + y2) * .5, 0);
		}

		userData = usd;
	};

	CurveDataAlum(QGEO::geoPoint3d ts, QGEO::geoPoint3d te, DWORD_PTR usd = NULL)
	{
		ps = ts;
		pe = te;

		userData = usd;

		isArc = false;

		pm.x = .5 * (ps.x + pe.x);
		pm.y = .5 * (ps.y + pe.y);
		pm.z = .5 * (ps.z + pe.z);
	};

	CurveDataAlum(QGEO::geoPoint3d ts, QGEO::geoPoint3d tm, QGEO::geoPoint3d te, DWORD_PTR usd = NULL)
	{
		ps = ts;
		pm = tm;
		pe = te;

		userData = usd;
		isArc = true;
	};

	void operator= (const CurveDataAlum& cv2)
	{
		ps = cv2.ps;
		pe = cv2.pe;
		pm = cv2.pm;

		userData = cv2.userData;
		isArc = cv2.isArc;
	};

	void Reverse()
	{
		CurveDataAlum cvn = CurveDataAlum(pe, ps, userData);
		*this = cvn;
	}
};

struct ARCHIOUT_LOOPDATA_ALUM
{
	std::vector<QGEO::geoPoint3d> m_pts;
	std::vector<CurveDataAlum> m_curves;
	QGEO::geoVector3d m_vecNorm; // 环 顺时针 还是 逆时针
	int iParentLoop;

	ARCHIOUT_LOOPDATA_ALUM()
	{};

	ARCHIOUT_LOOPDATA_ALUM(std::vector<CurveDataAlum> cvs)
	{
		m_curves.clear();
		m_curves = cvs;

		m_pts.clear();
		size_t len = m_curves.size();
		for (int i = 0; i < len; i++)
		{
			m_pts.push_back(m_curves[i].ps);
			if (i == len - 1)
			{
				m_pts.push_back(m_curves[i].pe);
			}
		}
	};

	ARCHIOUT_LOOPDATA_ALUM(std::vector<QGEO::geoPoint3d> pts)
	{
		// 封闭化处理
		size_t len = pts.size();
		if (len > 0)
		{
			QGEO::geoPoint3d p0 = pts[0];
			if (p0.DistanceTo(pts[len - 1]) > 1.0)
			{
				pts.push_back(p0);
			}
		}


		m_pts.clear();
		m_pts = pts;

		m_curves.clear();



		len = pts.size();
		if (len == 0)
		{
			return;
		}

		for (int i = 0; i < len - 1; i++)
		{
			int idnex = (i + 1) % len;
			CurveDataAlum cv(pts[i], pts[idnex]);
			m_curves.push_back(cv);
		}
	};


	DWORD_PTR userData;
};
typedef std::vector<ARCHIOUT_LOOPDATA_ALUM> loopArray;
typedef std::vector<CurveDataAlum> curveArray;

struct loopWithHole
{
	ARCHIOUT_LOOPDATA_ALUM lpOut;
	loopArray lpsHole;


	loopWithHole()
	{
		lpsHole.clear();
	};

	loopWithHole(ARCHIOUT_LOOPDATA_ALUM lp)
	{
		lpsHole.clear();
		lpOut = ARCHIOUT_LOOPDATA_ALUM(lp.m_pts);
		lpOut.userData = lp.userData;
	};

	loopWithHole(ARCHIOUT_LOOPDATA_ALUM lp, loopArray lpsHole2)
	{
		lpsHole.clear();
		lpOut = ARCHIOUT_LOOPDATA_ALUM(lp.m_pts);
		lpOut.userData = lp.userData;

		lpsHole.insert(lpsHole.end(), lpsHole2.begin(), lpsHole2.end());
	};

	void operator=(loopWithHole& ppl)
	{
		lpOut = ppl.lpOut;

		lpsHole.clear();
		lpsHole.insert(lpsHole.end(), ppl.lpsHole.begin(), ppl.lpsHole.end());
	};
};

typedef std::vector<loopWithHole> loopWithHoleArray;


struct CurveDataWithVecOut
{
	CurveDataAlum cv;
	QGEO::geoVector3d vecOut;

	/////// 新加的参数，适应跃层需要////////////////
	double topOffset;
	double botOffset;
	/////////////////////////////////////////////



	CurveDataWithVecOut()
	{
		topOffset = 0;
		botOffset = 0;
	};

	CurveDataWithVecOut(CurveDataAlum cv2, QGEO::geoVector3d vecOut2)
	{
		topOffset = 0;
		botOffset = 0;

		cv = cv2;
		vecOut = vecOut2;
	};
};

typedef std::vector<CurveDataWithVecOut> curveWithVecOutArray;


class CSDATADEF_EXPORT GraMath
{


private:
	static void GetNerestPjOfTwoPjs(CurveDataAlum& cv1, CurveDataAlum& cv2, geoPoint3d pj1, geoPoint3d pj2, QGEO::geoPoint3d& pj);
	static qCadDbObjectId GetLayerIDMust(const CString& layerName, COLORREF layerColor = RGB(255, 0, 0), short depth = 0, qCadDb::LineWeight lweight = qCadDb::qLnWt200);
	static bool IsArcEndPointsSameAsStartEndAngle(QGEO::geoPoint3d ptStart, QGEO::geoPoint3d ptMid, QGEO::geoPoint3d ptEnd);
	static bool GetOffsetDirection(QGEO::geoVector3d normal, bool bInside);
	static void AddCurveToFindLoop(ISearchLoop& searchLoop, CurveDataAlum cv);
	static void AddArcToFindLoop(ISearchLoop& searchLoop, geoArc arc, DWORD_PTR userData);
public:

	/*!
	 * @brief  构造函数
	*/
	GraMath();

	/*!
	 * @brief  析构函数
	*/
	~GraMath();

	/*!
	 * @brief vec 相等
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static void vecEq(std::vector<T>& v1, std::vector<T> v2)
	{
		v1.clear();
		v1.insert(v1.end(), v2.begin(), v2.end());
	};


	/*!
	 * @brief std::vector 反转
	 * @param v1
	*/
	template <typename T>
	static void vecReverse(std::vector<T>& v1)
	{
		std::reverse(v1.begin(), v1.end());
	};

	/*!
	 * @brief 排序的方法 去重
	 * @param v1
	*/
	template <typename T>
	static void vecUniqueSort(std::vector<T>& v1)
	{
		std::sort(v1.begin(), v1.end());
		v1.erase(std::unique(v1.begin(), v1.end()), v1.end());
	};


	/*!
	 * @brief 排序
	 * @param v1
	*/
	template <typename T>
	static void vecSort(std::vector<T>& v1)
	{
		std::sort(v1.begin(), v1.end());
	};

	/*!
 * @brief vecAppend
 * @param v1
 * @param v2
*/
	template <typename T>
	static void vecAppend(std::vector<T>& v1, std::vector<T> v2)
	{
		v1.insert(v1.end(), v2.begin(), v2.end());
	};


	/*!
	 * @brief vecHas
	 * @param v1
	 * @param da
	*/
	template <typename T>
	static bool vecHas(std::vector<T>& v1, T da)
	{
		auto iter = std::find(v1.begin(), v1.end(), da);
		if (iter != v1.end())
		{
			return true;
		}

		return false;
	};


	/*!
	 * @brief vecSetAt
	 * @param v1
	 * @param id
	 * @param da
	*/
	template <typename T>
	static void vecSetAt(std::vector<T>& v1, int id, T da)
	{
		int len = v1.size();
		for (int i = 0; i < len; i++)
		{
			if (id == i)
			{
				v1[i] = da;
				break;
			}
		}
	};


	/*!
	 * @brief vecUnion
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static std::vector<T> vecUnion(std::vector<T> v1, std::vector<T> v2)
	{
		std::set<T> s1(v1.begin(), v1.end());
		std::set<T> s2(v2.begin(), v2.end());

		std::set<T> s3 = lset_union(s1, s2);
		return lset2vec(s3);
	};


	/*!
	 * @brief vecDiff
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static std::vector<T> vecDiff(std::vector<T> v1, std::vector<T> v2)
	{
		std::set<T> s1(v1.begin(), v1.end());
		std::set<T> s2(v2.begin(), v2.end());

		std::set<T> s3 = lset_difference(s1, s2);
		return lset2vec(s3);
	};


	/*!
	 * @brief vecInter
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static std::vector<T> vecInter(std::vector<T> v1, std::vector<T> v2)
	{
		std::set<T> s1(v1.begin(), v1.end());
		std::set<T> s2(v2.begin(), v2.end());

		std::set<T> s3 = lset_intersection(s1, s2);
		return lset2vec(s3);
	};


	/*!
	 * @brief lsetHas
	 * @param v1
	 * @param da
	*/
	template <typename T>
	static bool lsetHas(std::set<T>& v1, T da)
	{
		auto iter = std::find(v1.begin(), v1.end(), da);
		if (iter != v1.end())
		{
			return true;
		}

		return false;
	};


	/*!
	 * @brief lsetErase
	 * @param v1
	 * @param da
	*/
	template <typename T>
	static bool lsetErase(std::set<T>& v1, T da)
	{
		auto iter = std::find(v1.begin(), v1.end(), da);
		if (iter != v1.end())
		{
			v1.erase(iter);
			return true;
		}

		return false;

	};


	/*!
	 * @brief lset_union
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static std::set<T> lset_union(std::set<T> v1, std::set<T> v2)
	{
		std::set<T> v3;
		v3.clear();
		set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), insert_iterator<std::set<int>>(v3, v3.begin()));
		return v3;
	};


	/*!
	 * @brief lset_intersection
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static std::set<T> lset_intersection(std::set<T> v1, std::set<T> v2)
	{
		std::set<T> v3;
		v3.clear();
		set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), insert_iterator<std::set<int>>(v3, v3.begin()));
		return v3;
	};


	/*!
	 * @brief lset_difference
	 * @param v1
	 * @param v2
	*/
	template <typename T>
	static std::set<T> lset_difference(std::set<T> v1, std::set<T> v2)
	{
		std::set<T> v3;
		v3.clear();
		set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), insert_iterator<std::set<int>>(v3, v3.begin()));
		return v3;
	};


	/*!
	 * @brief lset2vec
	 * @param s1
	*/
	template <typename T>
	static std::vector<T> lset2vec(std::set<T> s1)
	{
		std::vector<T> v1;
		v1.assign(s1.begin(), s1.end());

		return v1;
	};


	/*!
	 * @brief vec2lset
	 * @param v1
	*/
	template <typename T>
	static std::set<T> vec2lset(std::vector<T> v1)
	{
		std::set<T> s1(v1.begin(), v1.end());

		return s1;
	};



	/*!
	 * @brief mapHas
	 * @param mapDatas
	 * @param k1
	 * @param da
	*/
	template<typename T_key, typename T_da>
	static bool mapHas(std::map<T_key, T_da> mapDatas, T_key k1, T_da& da)
	{
		auto iter = mapDatas.find(k1);
		if (iter != mapDatas.end())
		{
			da = iter->second;
			return true;
		}

		return false;
	};



	/*!
	 * @brief mapInsert
	 * @param mapDatas
	 * @param k1
	 * @param da
	*/
	template<typename T_key, typename T_da>
	static void mapInsert(std::map<T_key, T_da>& mapDatas, T_key k1, T_da da)
	{
		T_da da0;
		if (mapHas(mapDatas, k1, da0))
		{
			mapDatas[k1] = da;
		}
		else
		{
			mapDatas.insert({ k1, da });
		}
	};




	/*!
	 * @brief Reverse
	 * @param t1
	 * @param t2
	*/
	template <typename T>
	static void Reverse(T& t1, T& t2)
	{
		std::swap(t1, t2);
	};


	/*!
	 * @brief vecRemoveAt
	 * @param v1
	 * @param id
	*/
	template <typename T>
	static void vecRemoveAt(std::vector<T>& v1, int id)
	{
		v1.erase(v1.begin() + id);
	};


	/*!
	 * @brief FromVector
	 * @param v1
	*/
	template <typename T>
	static void** FromVector(std::vector<T> v1)
	{
		T** pSideC = new T * [v1.size()];
		for (int i = 0; i < v1.size(); i++)
		{
			T* pC1 = new T();
			*pC1 = v1[i];
			pSideC[i] = pC1;
		}

		return (void**)&pSideC[0];
	};


	/*!
	 * @brief FromVector
	 * @param v1
	*/
	template <typename T>
	static void** FromVector(T v1)
	{
		T** pSideC = new T * [1];
		for (int i = 0; i < 1; i++)
		{
			T* pC1 = new T();
			*pC1 = v1;
			pSideC[i] = pC1;
		}

		return (void**)&pSideC[0];
	};

	static void ClearVectorMem(int len, void** p)
	{
		if (p)
		{
			for (int i = 0; i < len; i++)
			{
				if (p[i])
				{
					delete p[i];
					p[i] = NULL;
				}
			}

			delete p;
		}
	};


	/*!
	 * @brief ClearArrayP
	 * @param ars
	*/
	template <typename T>
	static void ClearArrayP(std::vector<T>& ars)
	{
		for (int i = 0; i < ars.size(); i++)
		{
			if (ars[i])
			{
				delete ars[i];
				ars[i] = NULL;
			}
		}

		ars.clear();
	};



	static QGEO::geoVector3d vec2d(QGEO::geoPoint3d ptFrom, QGEO::geoPoint3d ptTo);
	static QGEO::geoVector3d vec2d(CurveDataAlum cv);
	static QGEO::geoVector3d vec3d(QGEO::geoPoint3d ptFrom, QGEO::geoPoint3d ptTo);

	static std::vector<qCadDbObjectId> DrawPointsInGroup(geoPoint3dArray pts, qCadDbObjectId layerCurve);
	static qCadDbObjectId DrawPointInGroup(geoPoint3d pc, double r, qCadDbObjectId layerCurve);

	/*!
	 * @brief GetLayerID
	 * @param CString&
	 * @param layerColor
	*/
	static qCadDbObjectId GetLayerID(const CString& layerName, COLORREF layerColor = RGB(255, 0, 0), short depth = 0, qCadDb::LineWeight lweight = qCadDb::qLnWt200);

	static CurveDataAlum OffsetCurve(CurveDataAlum cv, QGEO::geoVector3d vecOff, double offset);
	static CurveDataAlum OffsetCurve(QGEO::geoPoint3d ptStart, QGEO::geoPoint3d ptMid, QGEO::geoPoint3d ptEnd, double offsetOut);

	static void GetTwoSideCurve(CurveDataAlum cvMid, double wid, CurveDataAlum& cvL, CurveDataAlum& cvR);
	static void GetTwoSideCurve(const geoPoint3d& ps, const geoPoint3d& pe, double wid, CurveDataAlum& cvL, CurveDataAlum& cvR);

	static ARCHIOUT_LOOPDATA_ALUM GetRectLoop(CurveDataAlum cvMid, double wid);
	static ARCHIOUT_LOOPDATA_ALUM GetRectLoop(CurveDataAlum cvL, CurveDataAlum cvR);
	static ARCHIOUT_LOOPDATA_ALUM GetRectLoop(const geoPoint3d& ps, const geoPoint3d& pe, double wid);
	static QGEO::geoVector3d vec2dLeft(CurveDataAlum cv);
	static QGEO::geoVector3d vec2dLeft(QGEO::geoPoint3d ms, QGEO::geoPoint3d me);
	static QGEO::geoVector3d vec2dRight(CurveDataAlum cv);
	static QGEO::geoVector3d vec2dRight(QGEO::geoPoint3d ms, QGEO::geoPoint3d me);
	static QGEO::geoPoint3d midPoint(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe);
	static void SetPtsZ(QGEO::geoPoint3dArray& pts, double elev);
	static void SetLoopElev(ARCHIOUT_LOOPDATA_ALUM& lp, double elev);
	static void SetLoopElev(loopArray& lps, double elev);
	static void SetLoopElevZero(ARCHIOUT_LOOPDATA_ALUM& lp);
	static QGEO::geoPoint3d PerToLine(QGEO::geoPoint3d pt, QGEO::geoPoint3d p1, QGEO::geoPoint3d p2);
	static QGEO::geoPoint3d PerToLine(QGEO::geoPoint3d pt, CurveDataAlum cv);
	static bool IsDVecEq(QGEO::geoVector3d v1, QGEO::geoVector3d v2, double tol = 0.001);
	static bool IsPointAtLeft(CurveDataAlum cvComp, geoPoint3d pt);
	static std::vector<QGEO::geoPoint3d> GetRectLoopPts(QGEO::geoPoint3d pc, double dSecL, double dSecW, double elev);
	static std::vector<qCadDbObjectId> DrawCurvesInGroup(curveArray cvs, qCadDbObjectId layerCurve, bool bBlock = false);
	static std::vector<qCadDbObjectId> DrawCurvesInGroupMust(curveArray cvs, qCadDbObjectId layerCurve, bool bBlock = false);
	static void ArcDataFrom3P(QGEO::geoPoint3d startPoint, QGEO::geoPoint3d pnt, QGEO::geoPoint3d endPoint, QGEO::geoPoint3d& m_center, double& m_radius, double& m_angleStart, double& m_angleEnd);
	static double AngleOfDVec3d(QGEO::geoVector3d vec);
	static double Distance3d(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe);
	static double Distance2d(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe);
	static double Distance2d(CurveDataAlum cv);
	static bool GetInterPointExt(CurveDataAlum cv1, CurveDataAlum cv2, QGEO::geoPoint3d& pj);
	static void CurveExtendBothSide(CurveDataAlum& cv, int disExt);
	static bool GetInterPointLine2LineExt(CurveDataAlum cvLin1, CurveDataAlum cvLin2, QGEO::geoPoint3d& pj);
	static void ArcToPtListEqDis(QGEO::geoPoint3d pc, float r, float as, float ae, std::vector<QGEO::geoPoint3d>& pts, int EqDis);
	static void MakePtsClosed(std::vector<QGEO::geoPoint3d>& pts, double dte = 2.);
	static INT64 GetInt(double data, int dte = 1);
	static QGEO::geoPoint3d FromSumOf(QGEO::geoPoint3d ps, QGEO::geoVector3d vec);
	static QGEO::geoPoint3d FromSumOf(QGEO::geoPoint3d ps, double ang, double dis);
	static std::vector<qCadDbObjectId> DrawLoopsInGroup(loopArray lps, qCadDbObjectId layerLoop, int offdisOut, bool bBlock = false);
	static qCadDbObjectId DrawPline(QGEO::geoPoint3dArray pts, qCadDbObjectId layerCurve);
	static void GetPtsExtend(QGEO::geoPoint3dArray pts, QGEO::geoPoint3d& pmin, QGEO::geoPoint3d& pmax);
	static CurveDataAlum curveFromArcData(geoPoint3d pc, double angs, double ange, double rad);
	static void SetCurveZ(CurveDataAlum& cv, double z);
	static bool IsDPointEq(QGEO::geoPoint3d p1, QGEO::geoPoint3d p2, double dte = 1.0);
	static bool IsDPointEq2d(QGEO::geoPoint3d p1, QGEO::geoPoint3d p2, double dte = 1.0);
	static QGEO::geoVector3d normalPolyGon(std::vector<QGEO::geoPoint3d> pts);
	static double GetPolygonDirection(QGEO::geoPoint3dArray& vtsIn);
	static ARCHIOUT_LOOPDATA_ALUM OffsetLoop(ARCHIOUT_LOOPDATA_ALUM lp, double wid, bool bInside);
	static BOOL GetCurveTrimdByStartCurve(CurveDataAlum& cv, CurveDataAlum cvPre);
	static void CurveTrimedByTwoSide(std::vector<CurveDataAlum>& curves, CurveDataAlum& cv, CurveDataAlum cvPre, CurveDataAlum cvNex, bool bLastCurve = false, bool bJoinPre = true, bool bJoinNex = true);
	static BOOL GetCurveTrimdByEndCurve(CurveDataAlum& cv, CurveDataAlum cvNex);
	static void OffsetLine(double offset, CurveDataAlum cv, CurveDataAlum& cvOff, bool bRight, bool bArcOutSide = true);
	static QGEO::geoVector3d RotateXY(QGEO::geoVector3d vec, double huAngle);
	static double Angle2d(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe);

	static void LineSplitedByLoopsVert(CurveDataAlum cv, std::vector<ARCHIOUT_LOOPDATA_ALUM> lps, std::vector<CurveDataAlum>& cvsOutLoop, std::vector<CurveDataAlum>& cvsInLoop, std::vector<CurveDataAlum>& cvsOnLoop);
	static curveArray To3d(curveArray cvs, geoMatrix3d mat);
	static CurveDataAlum To3d(CurveDataAlum cv, geoMatrix3d mat);
	static QGEO::geoPoint3d To3d(QGEO::geoPoint3d pt, geoMatrix3d mat);
	static void LineSplitedByLoops(CurveDataAlum cv, std::vector<ARCHIOUT_LOOPDATA_ALUM> lps, std::vector<CurveDataAlum>& cvsOutLoop, std::vector<CurveDataAlum>& cvsInLoop, std::vector<CurveDataAlum>& cvsOnLoop);
	static bool IsPointOnLine2d(QGEO::geoPoint3d pt, CurveDataAlum cv, double dte = 1.0);
	static bool IsPointOnLine2d(QGEO::geoPoint3d ps2, QGEO::geoPoint3d pe2, QGEO::geoPoint3d pt, double de = 1.0);
	static bool IsPointPerInLine(QGEO::geoPoint3d ps2, QGEO::geoPoint3d pe2, QGEO::geoPoint3d pt, double de = 1.0);
	static bool IsPointInLine(QGEO::geoPoint3d& pt, QGEO::geoPoint3d ps2, QGEO::geoPoint3d pe2, double de = 1.0);
	static double disToLine2d(QGEO::geoPoint3d pt, QGEO::geoPoint3d p1, QGEO::geoPoint3d p2);
	static std::vector<QGEO::geoPoint3d> OffsetPoints(std::vector<QGEO::geoPoint3d> pts, double wid, bool bInside);
	static void OffsetPoly(std::vector<double> widList, std::vector<QGEO::geoPoint3d> pts, std::vector<QGEO::geoPoint3d>& ptsNew, bool bInside);
	static void CurvesToPts(std::vector<CurveDataAlum> cvs, std::vector<QGEO::geoPoint3d>& pts);
	static void ValidPoint(QGEO::geoPoint3d& pt);
	static void ValidData(double& da);
	static QGEO::geoVector3d vec2dInner(CurveDataAlum cvSide, ARCHIOUT_LOOPDATA_ALUM loop, double dteInner = 10);
	static BOOL IsPointInPolyGon2d(QGEO::geoPoint3d pt, ARCHIOUT_LOOPDATA_ALUM lp);
	static BOOL IsPointInPolyGon2d(std::vector<QGEO::geoPoint3d> pts3, QGEO::geoPoint3d pt3d);
	static UINT PointInsidePolygon(std::vector<QGEO::geoPoint3d> points, const QGEO::geoPoint3d& pt, double EPSILON = 0.1);
	static bool IsPointInOneOfLoops(QGEO::geoPoint3d pm, std::vector<ARCHIOUT_LOOPDATA_ALUM> bounds, int offsetOutdis);
	static bool IsFaceHorizon(ARCHIOUT_LOOPDATA_ALUM lp);
	static bool IsFaceHorizon(const std::vector<QGEO::geoPoint3d>& vertexs);
	static void SetCurvesZ(curveArray& cvs, double z);
	static bool IsTwoLineParrallel(CurveDataAlum cv1, CurveDataAlum cv2, double dte = 0.01);
	static BOOL IsTwoLineParrallel(QGEO::geoPoint3d p1, QGEO::geoPoint3d p2, QGEO::geoPoint3d t1, QGEO::geoPoint3d t2, double dte = 0.01);
	static std::vector<CurveDataAlum> LineSplitedByCurves(CurveDataAlum cv, std::vector<CurveDataAlum> cvsCut);
	static ARCHIOUT_LOOPDATA_ALUM GetRectLoopVert(CurveDataAlum cv, double zmin, double zmax);
	static ARCHIOUT_LOOPDATA_ALUM GetRectLoopVert(const QGEO::geoPoint3d& ps, const QGEO::geoPoint3d& pe, double zmin, double zmax);
	static void SplitCurves(std::vector<CurveDataAlum> curves, std::vector<CurveDataAlum>& curvesResult, double tol_point_eq);
	static void ArcToPtList(QGEO::geoPoint3d pc, float r, float as, float ae, std::vector<QGEO::geoPoint3d>& pts, int nSplit);
	static CurveDataAlum To2d(CurveDataAlum cv, geoMatrix3d mat);
	static QGEO::geoPoint3d To2d(QGEO::geoPoint3d pt, geoMatrix3d mat);
	static loopArray To2d(loopArray lps, geoMatrix3d mat);
	static ARCHIOUT_LOOPDATA_ALUM To2d(ARCHIOUT_LOOPDATA_ALUM lp, geoMatrix3d mat);
	static geoMatrix3d matrixFrom(std::vector<QGEO::geoPoint3d> points);
	static geoMatrix3d matrixFrom(ARCHIOUT_LOOPDATA_ALUM lp);
	static BOOL IsFileExist(CString strFileName);
	static CString GetFileName(const CString& filePathName);
	static bool GetInterPointLine2ArcExt(CurveDataAlum cvLine, CurveDataAlum cvArc, QGEO::geoPoint3d& pj);

	static bool IsDVecParallel(QGEO::geoVector3d v1, QGEO::geoVector3d v2, double tol = 0.001);
	/*!
	 * @brief GetCentroid
	 * @param std::vector<QGEO::geoPoint3d>&
	*/
	static QGEO::geoPoint3d GetCentroid(const std::vector<QGEO::geoPoint3d>& vertexs);

	/*!
	 * @brief GetCentroid2d
	 * @param std::vector<QGEO::geoPoint3d>&
	*/
	static QGEO::geoPoint3d GetCentroid2d(const std::vector<QGEO::geoPoint3d>& vertexs);
	static BOOL IsStrHasStr(CString cs, CString csKey);
	static BOOL IsStrHasStr(CString cs, CString cs1, int& idhead, int& idtail);

	static void SetDrawInDrag(bool bDrag);
	static bool GetDrawInDrag();
};
