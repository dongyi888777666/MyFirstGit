#include "stdafx.h"
#include "GraMath.h"
#include "qCadDbCircle.h"
#include "qCadDbArc.h"
#include "dtlgeline2d.h"
#include "geoCircle.h"

bool g_bInDrag = false;

GraMath::GraMath()
{
}

GraMath::~GraMath()
{
}

QGEO::geoVector3d GraMath::vec2d(QGEO::geoPoint3d ptFrom, QGEO::geoPoint3d ptTo)
{
	QGEO::geoPoint3d ps2d(ptFrom.x, ptFrom.y, 0);
	QGEO::geoPoint3d pe2d(ptTo.x, ptTo.y, 0);


	return (pe2d - ps2d).Normalize();
}

QGEO::geoVector3d GraMath::vec2d(CurveDataAlum cv)
{
	return vec2d(cv.ps, cv.pe);
}

QGEO::geoVector3d GraMath::vec3d(QGEO::geoPoint3d ptFrom, QGEO::geoPoint3d ptTo)
{
	return (ptTo - ptFrom).Normalize();
}

std::vector<qCadDbObjectId> GraMath::DrawPointsInGroup(geoPoint3dArray pts, qCadDbObjectId layerCurve)
{
	std::vector<qCadDbObjectId> objids;

	int idd = 0;
	for (auto& pc : pts)
	{
		double r = 5 * (idd + 1);
		qCadDbObjectId idArc = DrawPointInGroup(pc, r, layerCurve);
		objids.push_back(idArc);
		idd++;
	}

	return objids;
}

qCadDbObjectId GraMath::DrawPointInGroup(geoPoint3d pc, double r, qCadDbObjectId layerCurve)
{
	qCadDbObjectId idArc;
//#ifdef DY_DEBUG_SHOW
	qCadDbCirclePtr pArc = qCadDbCircle::createObject();

	pArc->setCenter(pc);
	pArc->setRadius(r);
	pArc->setLayer(layerCurve);

	idArc = qCadInterfaceMgr::CreateNewEntity(pArc, layerCurve);
//#endif
	return idArc;
}

qCadDbObjectId GraMath::GetLayerIDMust(const CString& layerName, COLORREF layerColor, short depth, qCadDb::LineWeight lweight)
{
	qCadDbLayerTablePtr pTableLayer = qCadInterfaceMgr::getWorkingLayerTable();
	if (pTableLayer->has(layerName))
	{
		return pTableLayer->getAt(layerName);
	}
	qCadDbLayerTableRecordPtr  pNewLayer = qCadDbLayerTableRecord::createObject();
	qCadDbObjectId idNewLayer = pTableLayer->add(pNewLayer);
	pNewLayer->setName(layerName);
	pNewLayer->setColor(layerColor);
	///pNewLayer->setLinetypeObjectId(idLinetype);
	pNewLayer->setLineWeight(lweight);
	//pNewLayer->setLINETYPESCALE(100);//此处设置不起作用
	pTableLayer->close();
	return idNewLayer;
}

qCadDbObjectId GraMath::GetLayerID(const CString& layerName, COLORREF layerColor, short depth, qCadDb::LineWeight lweight)
{
//#ifdef DY_DEBUG_SHOW
	return GetLayerIDMust(layerName, layerColor, depth, lweight);
//#endif

	return qCadDbObjectId();
}

CurveDataAlum GraMath::OffsetCurve(CurveDataAlum cv, QGEO::geoVector3d vecOff, double offset)
{
	CurveDataAlum cvOff(cv.ps + vecOff * offset, cv.pe + vecOff * offset);
	return cvOff;
}

void GraMath::GetTwoSideCurve(const geoPoint3d& ps, const geoPoint3d& pe, double wid, CurveDataAlum& cvL, CurveDataAlum& cvR)
{
	CurveDataAlum cvMid(ps, pe);
	GetTwoSideCurve(cvMid, wid, cvL, cvR);
}


void GraMath::GetTwoSideCurve(CurveDataAlum cvMid, double wid, CurveDataAlum& cvL, CurveDataAlum& cvR)
{
	geoVector3d vLeft = GraMath::vec2dLeft(cvMid);
	cvL = GraMath::OffsetCurve(cvMid, vLeft, wid * .5);
	cvR = GraMath::OffsetCurve(cvMid, -vLeft, wid * .5);
}


ARCHIOUT_LOOPDATA_ALUM GraMath::GetRectLoop(const geoPoint3d& ps, const geoPoint3d& pe, double wid)
{
	CurveDataAlum cvMid(ps, pe);
	return GetRectLoop(cvMid, wid);
}

ARCHIOUT_LOOPDATA_ALUM GraMath::GetRectLoop(CurveDataAlum cvMid, double wid)
{
	CurveDataAlum cvL, cvR;
	GetTwoSideCurve(cvMid, wid, cvL, cvR);

	return GetRectLoop(cvL, cvR);
}


ARCHIOUT_LOOPDATA_ALUM GraMath::GetRectLoop(CurveDataAlum cvL, CurveDataAlum cvR)
{
	return ARCHIOUT_LOOPDATA_ALUM({ cvL.ps, cvL.pe, cvR.pe, cvR.ps, cvL.ps });
}

QGEO::geoVector3d GraMath::vec2dLeft(CurveDataAlum cv)
{
	return vec2dLeft(cv.ps, cv.pe);
}

QGEO::geoVector3d GraMath::vec2dLeft(QGEO::geoPoint3d ms, QGEO::geoPoint3d me)
{
	QGEO::geoVector3d vecHor = (me - ms).Normalize();
	return vecHor.CrossProduct(QGEO::geoVector3d::UNIT_AZ);
}

QGEO::geoVector3d GraMath::vec2dRight(CurveDataAlum cv)
{
	return vec2dRight(cv.ps, cv.pe);
}

QGEO::geoVector3d GraMath::vec2dRight(QGEO::geoPoint3d ms, QGEO::geoPoint3d me)
{
	QGEO::geoVector3d vecHor = (me - ms).Normalize();
	return vecHor.CrossProduct(QGEO::geoVector3d::UNIT_Z);
}

QGEO::geoPoint3d GraMath::midPoint(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe)
{
	double xm = (ps.x + pe.x) / 2;
	double ym = (ps.y + pe.y) / 2;
	double zm = (ps.z + pe.z) / 2;

	return QGEO::geoPoint3d(xm, ym, zm);
}


void GraMath::SetPtsZ(QGEO::geoPoint3dArray& pts, double elev)
{
	std::vector<QGEO::geoPoint3d> ptsNew;
	for each (auto pt in pts)
	{
		pt.z = elev;
		ptsNew.push_back(pt);
	}

	vecEq(pts, ptsNew);
}


void GraMath::SetLoopElev(ARCHIOUT_LOOPDATA_ALUM& lp, double elev)
{
	std::vector<QGEO::geoPoint3d> ptsNew;

	for each (auto pt in lp.m_pts)
	{
		pt.z = elev;
		ptsNew.push_back(pt);
	}

	lp = ARCHIOUT_LOOPDATA_ALUM(ptsNew);
}

void GraMath::SetLoopElev(loopArray& lps, double elev)
{
	for (auto iter = lps.begin(); iter != lps.end(); ++iter)
	{
		SetLoopElev((*iter), elev);
	}
}

void GraMath::SetLoopElevZero(ARCHIOUT_LOOPDATA_ALUM& lp)
{
	std::vector<QGEO::geoPoint3d> ptsNew;

	for each (auto pt in lp.m_pts)
	{
		pt.z = 0;
		ptsNew.push_back(pt);
	}

	lp = ARCHIOUT_LOOPDATA_ALUM(ptsNew);
}

QGEO::geoPoint3d GraMath::PerToLine(QGEO::geoPoint3d pt, CurveDataAlum cv)
{
	return PerToLine(pt, cv.ps, cv.pe);
}

QGEO::geoPoint3d GraMath::PerToLine(QGEO::geoPoint3d pt, QGEO::geoPoint3d p1, QGEO::geoPoint3d p2)
{
	geoLine l1 = geoLine(p1, p2);

	QGEO::geoPoint3d per;
	l1.GetPerpPt(pt, per);

	per.z = p1.z;

	return per;
}

bool GraMath::IsDVecEq(QGEO::geoVector3d v1, QGEO::geoVector3d v2, double tol /*= 0.001*/)
{
	geoTolerance gtol = geoTolerance(tol, tol);

	QGEO::geoVector3d vv1(v1.x, v1.y, v1.z);
	QGEO::geoVector3d vv2(v2.x, v2.y, v2.z);

	return vv1.isEqualTo(vv2, gtol);
}


bool GraMath::IsPointAtLeft(CurveDataAlum cvComp, geoPoint3d pt)
{
	geoPoint3d per = PerToLine(pt, cvComp.ps, cvComp.pe);

	geoVector3d vToPt = (per - pt);
	if (vToPt.Length() < 1) return false;

	vToPt = GraMath::vec2d(per, pt);

	geoVector3d vLeft = GraMath::vec2dLeft(cvComp);

	if (GraMath::IsDVecEq(vLeft, vToPt)) return true;
	return false;
}


std::vector<QGEO::geoPoint3d> GraMath::GetRectLoopPts(QGEO::geoPoint3d pc, double dSecL, double dSecW, double elev)
{
	std::vector<QGEO::geoPoint3d> ptVerts;
	ptVerts.emplace_back(QGEO::geoPoint3d(pc.x + dSecL / 2, pc.y + dSecW / 2, elev));
	ptVerts.emplace_back(QGEO::geoPoint3d(pc.x - dSecL / 2, pc.y + dSecW / 2, elev));
	ptVerts.emplace_back(QGEO::geoPoint3d(pc.x - dSecL / 2, pc.y - dSecW / 2, elev));
	ptVerts.emplace_back(QGEO::geoPoint3d(pc.x + dSecL / 2, pc.y - dSecW / 2, elev));

	return ptVerts;
}

std::vector<qCadDbObjectId> GraMath::DrawCurvesInGroup(curveArray cvs, qCadDbObjectId layerCurve, bool bBlock)
{
	std::vector<qCadDbObjectId> objids;
	objids = DrawCurvesInGroupMust(cvs, layerCurve, bBlock);
	return objids;
}

std::vector<qCadDbObjectId> GraMath::DrawCurvesInGroupMust(curveArray cvs, qCadDbObjectId layerCurve, bool bBlock)
{
	std::vector<qCadDbObjectId> objids;

	if (bBlock)
	{
		//CreateBlock(ents);
	}
	else
	{
		for each (auto & cv in cvs)
		{
			if (cv.isArc)
			{
				qCadDbArcPtr pArc = qCadDbArc::createObject();

				QGEO::geoPoint3d pc;
				double r, as, ae;
				ArcDataFrom3P(cv.ps, cv.pm, cv.pe, pc, r, as, ae);


				pArc->setCenter(pc);
				pArc->setRadius(r);
				pArc->setStartAngle(as);
				pArc->setEndAngle(ae);
				pArc->setLayer(layerCurve);

				qCadDbObjectId idArc = qCadInterfaceMgr::CreateNewEntity(pArc, layerCurve);
				objids.push_back(idArc);
			}
			else
			{
				if (GraMath::Distance3d(cv.ps, cv.pe) > 2)
				{
					qCadDbLinePtr pLine = qCadDbLine::createObject();
					pLine->setStartPoint(cv.ps);
					pLine->setEndPoint(cv.pe);
					pLine->setLayer(layerCurve);
					pLine->setColorIndex(256);	//0--随块，256--随层		

					qCadDbObjectId idLine = qCadInterfaceMgr::CreateNewEntity(pLine, layerCurve);
					objids.push_back(idLine);
				}
			}
		}
	}
	return objids;
}


void GraMath::ArcDataFrom3P(QGEO::geoPoint3d startPoint, QGEO::geoPoint3d pnt, QGEO::geoPoint3d endPoint, QGEO::geoPoint3d& m_center, double& m_radius, double& m_angleStart, double& m_angleEnd)
{
	m_angleStart = 0;
	m_angleEnd = 0;
	m_radius = 1;

	if (startPoint.isEqualTo(endPoint, 1.0))
	{
		m_angleStart = 0;
		m_angleEnd = M_PI * 2;
		m_radius = startPoint.DistanceTo(pnt) / 2;
		m_center = (startPoint + pnt) / 2;

		m_center.z = (startPoint.z + endPoint.z) * .5;

		return;
	}

	QGEO::geoVector3d v11 = (startPoint - endPoint).Normalize();
	QGEO::geoVector3d v2 = (pnt - endPoint).Normalize();


	if ((v11 - v2).IsZeroLength() || (v11 + v2).IsZeroLength())
	{
		m_angleStart = 0;
		m_angleEnd = M_PI * 2;
		m_radius = startPoint.DistanceTo(pnt) / 2;
		m_center = (startPoint + pnt) / 2;
		m_center.z = (startPoint.z + endPoint.z) * .5;
		return;
	}
	double u1 = (pnt.x * pnt.x - startPoint.x * startPoint.x +
		pnt.y * pnt.y - startPoint.y * startPoint.y) / 2.0;
	double u2 = (endPoint.x * endPoint.x - startPoint.x * startPoint.x +
		endPoint.y * endPoint.y - startPoint.y * startPoint.y) / 2.0;
	double d11 = pnt.x - startPoint.x;
	double d12 = pnt.y - startPoint.y;
	double d21 = endPoint.x - startPoint.x;
	double d22 = endPoint.y - startPoint.y;
	double denom = (d11 * d22 - d21 * d12);
	m_center.x = (u1 * d22 - u2 * d12) / denom;
	m_center.y = (u2 * d11 - u1 * d21) / denom;
	m_center.z = (startPoint.z + endPoint.z) * .5;

	m_radius = m_center.DistanceTo(startPoint);
	m_angleStart = AngleOfDVec3d(startPoint - m_center); // 此處可能有問題
	m_angleEnd = AngleOfDVec3d(endPoint - m_center);
	double angleMid = AngleOfDVec3d(pnt - m_center);
	BOOL b1 = (angleMid < m_angleStart) && (angleMid > m_angleEnd);
	BOOL b01 = (angleMid > m_angleStart) && (angleMid > m_angleEnd);
	BOOL b02 = (angleMid < m_angleStart) && (angleMid < m_angleEnd);
	BOOL b2 = (b01 || b02) && (m_angleStart < m_angleEnd);
	if (b1 || b2) std::swap(m_angleStart, m_angleEnd);

	//return *this;
}

double GraMath::AngleOfDVec3d(QGEO::geoVector3d vec)
{
	return vec.AngleTo_uu(QGEO::geoVector3d::UNIT_X, QGEO::geoVector3d::UNIT_AZ);
}

double GraMath::Distance3d(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe)
{
	return ps.DistanceTo(pe);
}

double GraMath::Distance2d(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe)
{
	QGEO::geoPoint3d ps2d(ps.x, ps.y, 0);
	QGEO::geoPoint3d pe2d(pe.x, pe.y, 0);

	return ps2d.DistanceTo(pe2d);
}


double GraMath::Distance2d(CurveDataAlum cv)
{
	return Distance2d(cv.ps, cv.pe);
}

bool GraMath::GetInterPointExt(CurveDataAlum cv1, CurveDataAlum cv2, QGEO::geoPoint3d& pj)
{
	if (cv1.isArc)
	{
		geoArc arc1(cv1.ps, cv1.pm, cv1.pe);



		if (cv2.isArc)
		{
			return false;

			//return GetInterPointArc2ArcExt(cv1, cv2, pj); // 暂时关闭
		}
		else
		{
			return GetInterPointLine2ArcExt(cv2, cv1, pj);
		}
	}
	else
	{
		if (cv2.isArc)
		{
			return GetInterPointLine2ArcExt(cv1, cv2, pj);
		}
		else
		{
			return GetInterPointLine2LineExt(cv1, cv2, pj);
		}
	}

	return false;
}

bool GraMath::GetInterPointLine2ArcExt(CurveDataAlum cvLine, CurveDataAlum cvArc, QGEO::geoPoint3d& pj)
{
	gePoint2d ps1(cvArc.ps.x, cvArc.ps.y);
	gePoint2d pe1(cvArc.pe.x, cvArc.pe.y);

	gePoint2d ps2(cvLine.ps.x, cvLine.ps.y);
	gePoint2d pe2(cvLine.pe.x, cvLine.pe.y);


	CurveDataAlum cvLineExt = cvLine;
	GraMath::CurveExtendBothSide(cvLineExt, 2000);

	geoLine lin2(cvLineExt.ps, cvLineExt.pe);


	QGEO::geoPoint3d pc;
	double r, as, ae;
	ArcDataFrom3P(cvArc.ps, cvArc.pm, cvArc.pe, pc, r, as, ae);

	geoCircle cir(pc, r);


	if (abs(r - 5520) < 2)
	{
		int kkk = 0;
	}

	qAllocArray<QGEO::geoPoint3d> aryInsecPoints;
	geoTolerance eps; eps.setDisTol(1);
	int num = cir.insectionwithCurve(&lin2, aryInsecPoints, eps);

	if (num == 2)
	{
		GetNerestPjOfTwoPjs(cvArc, cvLine, aryInsecPoints[0], aryInsecPoints[1], pj);
		pj.z = cvArc.ps.z;
		return true;
	}
	else if (num == 1)
	{
		pj.Set(aryInsecPoints[0].x, aryInsecPoints[0].y, cvArc.ps.z);
		return true;
	}

	return false;
}


void GraMath::GetNerestPjOfTwoPjs(CurveDataAlum& cv1, CurveDataAlum& cv2, geoPoint3d pj1, geoPoint3d pj2, QGEO::geoPoint3d& pj)
{
	geoPoint3d ps1 = cv1.ps;
	geoPoint3d pe1 = cv1.pe;
	geoPoint3d ps2 = cv2.ps;
	geoPoint3d pe2 = cv2.pe;

	double dis1s2 = ps1.DistanceTo(ps2);
	double dis1e2 = ps1.DistanceTo(pe2);
	double die1s2 = pe1.DistanceTo(ps2);
	double die1e2 = pe1.DistanceTo(pe2);

	double mindis = min(dis1s2, min(dis1e2, min(die1s2, die1e2)));
	geoPoint3d pTest = ps1;
	if (mindis == dis1s2 || mindis == dis1e2)
	{
		pTest = ps1;
	}
	else
	{
		pTest = pe1;
	}

	double disj1 = pTest.DistanceTo(pj1);
	double disj2 = pTest.DistanceTo(pj2);

	if (disj1 < disj2)
	{
		pj = pj1;
	}
	else
	{
		pj = pj2;
	}
}

void GraMath::CurveExtendBothSide(CurveDataAlum& cv, int disExt)
{
	QGEO::geoVector3d vecHor = GraMath::vec2d(cv);
	QGEO::geoPoint3d ds = cv.ps;
	QGEO::geoPoint3d de = cv.pe;
	ds = ds - vecHor * disExt;
	de = de + vecHor * disExt;

	cv = CurveDataAlum(ds, de);
}

bool GraMath::GetInterPointLine2LineExt(CurveDataAlum cvLin1, CurveDataAlum cvLin2, QGEO::geoPoint3d& pj)
{
	CurveDataAlum cvLineExt1 = cvLin1;
	GraMath::CurveExtendBothSide(cvLineExt1, 2000);

	CurveDataAlum cvLineExt2 = cvLin2;
	GraMath::CurveExtendBothSide(cvLineExt2, 2000);

	gePoint2d ps1(cvLineExt1.ps.x, cvLineExt1.ps.y);
	gePoint2d pe1(cvLineExt1.pe.x, cvLineExt1.pe.y);

	gePoint2d ps2(cvLineExt2.ps.x, cvLineExt2.ps.y);
	gePoint2d pe2(cvLineExt2.pe.x, cvLineExt2.pe.y);



	geLine2d lin1(ps1, pe1);
	geLine2d lin2(ps2, pe2);

	gePoint2d pj1, pj2;
	int n1 = lin1.IntersectWithLine(lin2, pj1, pj2, FALSE, xTol());
	if (n1 == 1)
	{
		pj.Set(pj1.x, pj1.y, cvLineExt1.ps.z);
		return true;
	}

	return false;
}


void GraMath::ArcToPtListEqDis(QGEO::geoPoint3d pc, float r, float as, float ae, std::vector<QGEO::geoPoint3d>& pts, int EqDis)
{
	// 全部采用弧度制
	double EqAng = EqDis / r;
	int nSplit = GraMath::GetInt(2 * M_PI / EqAng);

	pts.clear();
	if (as > ae)
	{
		ae += 2 * M_PI;
	}

	for (size_t i = 0; i < nSplit + 2; i++)
	{
		double ang = as + i * EqAng;
		if (ang > ae) break;

		QGEO::geoVector3d vec(cos(ang) * r, sin(ang) * r, 0);
		QGEO::geoPoint3d pp = FromSumOf(pc, vec);
		pts.push_back(pp);
	}

	QGEO::geoVector3d vec(cos(ae) * r, sin(ae) * r, 0);
	QGEO::geoPoint3d pp = FromSumOf(pc, vec);
	pts.push_back(pp);

}

void GraMath::MakePtsClosed(std::vector<QGEO::geoPoint3d>& pts, double dte)
{
	int len = pts.size();
	if (len < 3) return;

	QGEO::geoPoint3d p0 = pts[0];
	if (p0.DistanceTo(pts[len - 1]) > dte)
	{
		pts.push_back(p0);
	}
}

INT64 GraMath::GetInt(double data, int dte /*= 1*/)
{
	double datam = fabs(data);


	double d = datam / (dte + 0.0);
	INT64 n = (INT64)d;
	if (d - n > 0.5) n++;

	INT64 iDat = n * dte;
	if (data < 0) iDat *= -1;

	return iDat;
}

QGEO::geoPoint3d GraMath::FromSumOf(QGEO::geoPoint3d ps, double ang, double dis)
{
	QGEO::geoVector3d vec(cos(ang) * dis, sin(ang) * dis, 0);
	QGEO::geoPoint3d pp = FromSumOf(ps, vec);

	return pp;
}

QGEO::geoPoint3d GraMath::FromSumOf(QGEO::geoPoint3d ps, QGEO::geoVector3d vec)
{
	QGEO::geoPoint3d pe = ps + vec;
	return pe;
}

std::vector<qCadDbObjectId> GraMath::DrawLoopsInGroup(loopArray lps, qCadDbObjectId layerLoop, int offdisOut, bool bBlock)
{
	std::vector<qCadDbObjectId> ids;
	for each (auto & lp in lps)
	{
		ids.push_back(DrawPline(lp.m_pts, layerLoop));
	}

	if (bBlock)
	{
		//ids.push_back(CreateBlock(entsBlock));
	}

	return ids;
}



qCadDbObjectId GraMath::DrawPline(QGEO::geoPoint3dArray pts, qCadDbObjectId layerCurve)
{
	std::vector<int> zs;
	bool bHor = true;
	for each (auto & pt in pts)
	{
		int zz = (int)pt.z;
		if (!vecHas(zs, zz))
		{
			zs.push_back(zz);

			if (zs.size() > 1)
			{
				bHor = false;
				break;
			}
		}
	}

	qCadDbObjectId lineId;
	//if (!bHor)
	//{
	//	// 不是水平面
	//	return DrawFaceByLayer(pts, layerCurve);
	//}


	QGEO::geoPoint3d pmin, pmax;
	GetPtsExtend(pts, pmin, pmax);

	if ((pmax.x - pmin.x) < 5 && (pmax.y - pmin.y) < 5) return lineId;


	int elev = 0;
	qCadDbPolylinePtr polyLine = qCadDbPolyline::createObject();
	for (int i = 0; i < pts.size(); i++)
	{
		QGEO::geoPoint3d pt = pts[i];
		geoPoint2d pto(pt.x, pt.y);
		elev = pt.z;
		polyLine->addVertexAt(i, pto);
	}

	polyLine->setClosed(false);
	polyLine->setElevation(elev);
	polyLine->setLineWeight(qCadDb::qLnWt200);




	if (layerCurve.isNull()) { lineId = qCadInterfaceMgr::CreateNewEntity(polyLine); }
	else { lineId = qCadInterfaceMgr::CreateNewEntity(polyLine, layerCurve); }

	return lineId;
}


void GraMath::GetPtsExtend(QGEO::geoPoint3dArray pts, QGEO::geoPoint3d& pmin, QGEO::geoPoint3d& pmax)
{
	if (pts.empty()) return;

	int len = pts.size();
	double minx, miny, maxx, maxy;
	for (int i = 0; i < len; i++)
	{
		if (i == 0)
		{
			minx = pts[i].x;
			miny = pts[i].y;
			maxx = pts[i].x;
			maxy = pts[i].y;
		}
		else
		{
			if (pts[i].x < minx) minx = pts[i].x;
			if (pts[i].y < miny) miny = pts[i].y;
			if (pts[i].x > maxx) maxx = pts[i].x;
			if (pts[i].y > maxy) maxy = pts[i].y;
		}
	}

	pmin.Set(minx, miny, pts[0].z);
	pmax.Set(maxx, maxy, pts[0].z);
}

CurveDataAlum GraMath::curveFromArcData(geoPoint3d pc, double angs, double ange, double rad)
{
	QGEO::geoPoint3d ps = GraMath::FromSumOf(pc, angs, rad);
	QGEO::geoPoint3d pe = GraMath::FromSumOf(pc, ange, rad);


	geoVector3d vecPath = GraMath::vec2d(ps, pe);
	QGEO::geoVector3d vecPer = vecPath.CrossProduct(QGEO::geoVector3d::UNIT_Z).GetNormal();

	double angDiff = (ange - angs);
	if (angDiff < 0) angDiff += (M_PI * 2);
	if (angDiff > M_PI * 2) angDiff -= (M_PI * 2);

	double angm = angs + angDiff * .5;
	QGEO::geoPoint3d pm = GraMath::FromSumOf(pc, angm, rad);

	return CurveDataAlum(ps, pm, pe);
}

void GraMath::SetCurveZ(CurveDataAlum& cv, double z)
{
	cv.ps.z = z;
	cv.pe.z = z;
	cv.pm.z = z;
}


bool GraMath::IsDPointEq(QGEO::geoPoint3d p1, QGEO::geoPoint3d p2, double dte)
{
	if (abs(p1.x - p2.x) > dte) return false;
	if (abs(p1.y - p2.y) > dte) return false;
	if (abs(p1.z - p2.z) > dte) return false;

	return true;
}

bool GraMath::IsDPointEq2d(QGEO::geoPoint3d p1, QGEO::geoPoint3d p2, double dte)
{
	if (abs(p1.x - p2.x) > dte) return false;
	if (abs(p1.y - p2.y) > dte) return false;

	return true;
}

QGEO::geoVector3d GraMath::normalPolyGon(std::vector<QGEO::geoPoint3d> pts)
{
	if (pts.size() < 3) return QGEO::geoVector3d(0, 0, 0);
	QGEO::geoPoint3d p1 = QGEO::geoPoint3d(pts[0].x, pts[0].y, pts[0].z);


	bool bFind2 = false;
	QGEO::geoPoint3d p2;
	int len = pts.size();

	int id2 = 1;
	for (int i = 1; i < len; i++)
	{
		p2 = pts[i];

		if (p1.DistanceTo(p2) > 2)
		{
			bFind2 = true;
			id2 = i;
			break;
		}
	}

	if (!bFind2) return QGEO::geoVector3d(0, 0, 0);

	QGEO::geoVector3d vec1 = (p2 - p1).Normalize();
	QGEO::geoVector3d vec2 = vec1;
	bool bFind3 = false;

	QGEO::geoPoint3d p3;
	for (int i = id2 + 1; i < len; i++)
	{
		p3 = pts[i];
		vec2 = (p3 - p2).Normalize();

		if (!GraMath::IsDVecParallel(vec1, vec2))
		{
			bFind3 = true;
			break;
		}
	}

	if (!bFind3) return QGEO::geoVector3d(0, 0, 0);

	QGEO::geoVector3d vecNorm = vec1.CrossProduct(vec2).Normalize();

	if (GraMath::IsDVecParallel(vecNorm, QGEO::geoVector3d::UNIT_Z))
	{
		// 判断多边形的逆时针还是顺时针
		double area = GetPolygonDirection(pts);

		if (area > 0) vecNorm.z = 1;
		else vecNorm.z = -1;
	}


	return vecNorm;
}


double GraMath::GetPolygonDirection(QGEO::geoPoint3dArray& vtsIn)
{
	if (vtsIn.empty()) return 0;


	bool bEq = vtsIn.back().isEqualTo(vtsIn.front(), 1.0);

	// if not closed then close it
	if (!bEq)
	{
		vtsIn.push_back(vtsIn[0]);
	}

	double area_polygon = 0;
	// computer the polyon area(twice of double value)
	for (UINT i = 1; i <= vtsIn.size() - 1; ++i)
	{
		area_polygon = vtsIn[i - 1].x * vtsIn[i].y - vtsIn[i].x * vtsIn[i - 1].y + area_polygon;
	}

	return  area_polygon;
}

bool GraMath::IsDVecParallel(QGEO::geoVector3d v1, QGEO::geoVector3d v2, double tol /*= 0.001*/)
{
	geoTolerance gtol = geoTolerance(tol, tol);

	if (v1.IsParallelTo(v2, gtol) == 0) return false;
	return true;
}

bool GraMath::GetOffsetDirection(QGEO::geoVector3d normal, bool bInside)
{
	if (normal.z > 0)
	{
		return !bInside;
	}
	else
	{
		return bInside;
	}
}

ARCHIOUT_LOOPDATA_ALUM GraMath::OffsetLoop(ARCHIOUT_LOOPDATA_ALUM lp, double wid, bool bInside)
{
	if (wid == 0) return lp;

	QGEO::geoVector3d normal = GraMath::normalPolyGon(lp.m_pts);
	bool bRight = GetOffsetDirection(normal, bInside);

	int ncv = lp.m_curves.size();

	curveArray cvsNew;
	for (int i = 0; i < ncv; i++)
	{
		int iPre = (i - 1 + ncv) % ncv;
		int iNex = (i + 1) % ncv;

		CurveDataAlum cv = lp.m_curves[i];
		CurveDataAlum cvPre = lp.m_curves[iPre];
		CurveDataAlum cvNex = lp.m_curves[iNex];

		CurveDataAlum cvOff, cvPreOff, cvNexOff;
		OffsetLine(wid, cv, cvOff, bRight);
		OffsetLine(wid, cvPre, cvPreOff, bRight);
		OffsetLine(wid, cvNex, cvNexOff, bRight);


		//if (ncv == 1)
		//{
		//	cvsNew.clear();
		//	cvsNew.push_back(cvOff);
		//	return;
		//}

		bool bLastCurve = (i == ncv - 1) ? true : false;
		bool bJoinPre = true, bJoinNex = true;



		CurveTrimedByTwoSide(cvsNew, cvOff, cvPreOff, cvNexOff, bLastCurve, bJoinPre, bJoinNex);
	}

	ARCHIOUT_LOOPDATA_ALUM lpNew = ARCHIOUT_LOOPDATA_ALUM(cvsNew);
	return lpNew;
}


BOOL GraMath::GetCurveTrimdByStartCurve(CurveDataAlum& cv, CurveDataAlum cvPre)
{
	int inum = 0;
	QGEO::geoPoint3d j1, j2;
	if (GraMath::GetInterPointExt(cv, cvPre, j1))
	{
		if (cv.isArc)
		{
			cv = CurveDataAlum(j1, cv.pm, cv.pe);
		}
		else
		{
			cv.ps = j1;
		}
		return TRUE;
	}

	return FALSE;
}



void GraMath::CurveTrimedByTwoSide(std::vector<CurveDataAlum>& curves, CurveDataAlum& cv, CurveDataAlum cvPre, CurveDataAlum cvNex, bool bLastCurve, bool bJoinPre, bool bJoinNex)
{
	if (bJoinPre)
	{
		if (!GetCurveTrimdByStartCurve(cv, cvPre))
		{
			// 和前一条线，是平行关系
			if (cvPre.pe.DistanceTo(cv.ps) > 1)
			{
				// 增加连接线段，保证轮廓封闭
				CurveDataAlum cvm(cvPre.pe, cv.ps, cv.userData);
				curves.push_back(cvm);
			}
		}
	}

	if (bJoinNex)
	{
		// 只考虑和前面的线条平行时，增加链接线段
		GetCurveTrimdByEndCurve(cv, cvNex);
	}

	curves.push_back(cv);

	//BOOL isParrallelPre = PublicFunc::IsTwoLineParrallel(cv.ps, cv.pe, cvPre.ps, cvPre.pe);
	//BOOL isParrallelNex = PublicFunc::IsTwoLineParrallel(cv.ps, cv.pe, cvNex.ps, cvNex.pe);

	//if (isParrallelPre)
	//{
	//	QGEO::geoPoint3d pj1, pj2, pj3;
	//	FindParrelLineLinkPoint(cvPre.ps, cvPre.pe, cv.ps, cv.pe, pj1, pj2, pj3);

	//	if (pj1.DistanceTo(pj2) > 1)
	//	{
	//		// 需要额外添加两段
	//		CurveDataAlum cv1(pj1, pj2);
	//		CurveDataAlum cv2(pj2, pj3);
	//		curves.push_back(cv1);
	//		curves.push_back(cv2);
	//	}
	//}
	//else
	//{
	//	if (bJoinPre)
	//	{
	//		QGEO::geoPoint3d pjs;
	//		PublicFunc::GetInterPoint(cv.ps, cv.pe, cvPre.ps, cvPre.pe, pjs);
	//		cv.ps = pjs;
	//	}
	//}

	//if (!isParrallelNex)
	//{
	//	if (bJoinNex)
	//	{
	//		QGEO::geoPoint3d pje;
	//		PublicFunc::GetInterPoint(cv.ps, cv.pe, cvNex.ps, cvNex.pe, pje);
	//		cv.pe = pje;
	//	}
	//}
	//
	//PublicFunc::GetMidPoint(cv.ps, cv.pe, cv.pm);
	//curves.push_back(cv);

	//if (bLastCurve&&isParrallelNex)
	//{
	//	QGEO::geoPoint3d pj1, pj2, pj3;
	//	FindParrelLineLinkPoint(cv.ps, cv.pe, cvNex.ps, cvNex.pe, pj1, pj2, pj3);
	//	if (pj1.DistanceTo(pj2) > 1)
	//	{
	//		// 需要额外添加两段
	//		CurveDataAlum cv1(pj1, pj2);
	//		CurveDataAlum cv2(pj2, pj3);
	//		curves.push_back(cv1);
	//		curves.push_back(cv2);
	//	}
	//}
}


BOOL GraMath::GetCurveTrimdByEndCurve(CurveDataAlum& cv, CurveDataAlum cvNex)
{
	int inum = 0;
	QGEO::geoPoint3d j1;
	if (GraMath::GetInterPointExt(cv, cvNex, j1))
	{
		if (cv.isArc)
		{
			cv = CurveDataAlum(cv.ps, cv.pm, j1);
		}
		else
		{
			cv.pe = j1;
		}
		return TRUE;
	}

	return FALSE;
}

void GraMath::OffsetLine(double offset, CurveDataAlum cv, CurveDataAlum& cvOff, bool bRight, bool bArcOutSide)
{
	if (fabs(offset) < 1.0)
	{
		cvOff = cv;
		return;
	}

	cvOff = cv;
	//cvOff.isArc = cv.isArc;
	if (cv.isArc == false)
	{
		QGEO::geoVector3d vec = (cv.pe - cv.ps).Normalize() * offset;

		if (bRight)
			vec = GraMath::RotateXY(vec, -M_PI / 2);
		else vec = GraMath::RotateXY(vec, M_PI / 2);

		cvOff.ps = GraMath::FromSumOf(cv.ps, vec);
		cvOff.pm = GraMath::FromSumOf(cv.pm, vec);
		cvOff.pe = GraMath::FromSumOf(cv.pe, vec);
	}
	else
	{
		if (bArcOutSide)
		{
			cvOff = GraMath::OffsetCurve(cv.ps, cv.pm, cv.pe, offset);
		}
		else
		{
			cvOff = GraMath::OffsetCurve(cv.ps, cv.pm, cv.pe, -offset);
		}
	}
}


CurveDataAlum GraMath::OffsetCurve(QGEO::geoPoint3d ptStart, QGEO::geoPoint3d ptMid, QGEO::geoPoint3d ptEnd, double offsetOut)
{
	QGEO::geoPoint3d pc;
	double r, as, ae;
	GraMath::ArcDataFrom3P(ptStart, ptMid, ptEnd, pc, r, as, ae);

	QGEO::geoVector3d vecs = (ptStart - pc).Normalize();
	QGEO::geoVector3d vece = (ptEnd - pc).Normalize();
	QGEO::geoVector3d vecm = (ptMid - pc).Normalize();

	double radNew = r + offsetOut;
	if (IsArcEndPointsSameAsStartEndAngle(ptStart, ptMid, ptEnd))
	{
		radNew = r - offsetOut;
	}

	QGEO::geoPoint3d ps = pc + vecs * radNew;
	QGEO::geoPoint3d pe = pc + vece * radNew;
	QGEO::geoPoint3d pm = pc + vecm * radNew;

	return CurveDataAlum(ps, pm, pe);
}

QGEO::geoVector3d GraMath::RotateXY(QGEO::geoVector3d vec, double huAngle)
{
	geoMatrix _tm;
	_tm.OfRotation(huAngle, QGEO::geoVector3d::UNIT_Z, geoPoint(0, 0, 0));

	QGEO::geoVector3d v = vec;
	v.TransformBy(_tm);
	return QGEO::geoVector3d(v.x, v.y, v.z);
}

double GraMath::Angle2d(QGEO::geoPoint3d ps, QGEO::geoPoint3d pe)
{
	geoVector3d vs = GraMath::vec2d(ps, pe);
	return AngleOfDVec3d(vs);
}


bool GraMath::IsArcEndPointsSameAsStartEndAngle(QGEO::geoPoint3d ptStart, QGEO::geoPoint3d ptMid, QGEO::geoPoint3d ptEnd)
{
	QGEO::geoPoint3d pc;
	double r, as, ae;
	GraMath::ArcDataFrom3P(ptStart, ptMid, ptEnd, pc, r, as, ae);

	double asReal = GraMath::Angle2d(pc, ptStart);
	double aeReal = GraMath::Angle2d(pc, ptEnd);



	double angsdiff = as - asReal;
	if (angsdiff > 2 * M_PI) angsdiff -= 2 * M_PI;


	double angediff = as - aeReal;
	if (angediff > 2 * M_PI) angediff -= 2 * M_PI;

	if (abs(angsdiff) < abs(angediff))
	{
		return true;
	}






	//double angdiff = as - asReal;
	//if (angdiff > 2 * M_PI) angdiff -= 2 * M_PI;

	//if (abs(angdiff) < 0.1)
	//{
	//	return true;
	//}

	return false;
}



void GraMath::LineSplitedByLoopsVert(CurveDataAlum cv, std::vector<ARCHIOUT_LOOPDATA_ALUM> lps,
	std::vector<CurveDataAlum>& cvsOutLoop,
	std::vector<CurveDataAlum>& cvsInLoop,
	std::vector<CurveDataAlum>& cvsOnLoop)
{
	if (lps.size() == 0) return;
	geoMatrix3d mat = GraMath::matrixFrom(lps[0].m_pts);

	CurveDataAlum cvPlan = To2d(cv, mat);
	loopArray lpsPlan = To2d(lps, mat);

	std::vector<CurveDataAlum> cvsOutLoopPlan;
	std::vector<CurveDataAlum> cvsInLoopPlan;
	std::vector<CurveDataAlum> cvsOnLoopPlan;
	GraMath::LineSplitedByLoops(cvPlan, lpsPlan, cvsOutLoopPlan, cvsInLoopPlan, cvsOnLoopPlan);

	cvsOutLoop.clear();
	GraMath::vecAppend(cvsOutLoop, To3d(cvsOutLoopPlan, mat));

	cvsInLoop.clear();
	GraMath::vecAppend(cvsInLoop, To3d(cvsInLoopPlan, mat));

	cvsOnLoop.clear();
	GraMath::vecAppend(cvsOnLoop, To3d(cvsOnLoopPlan, mat));
}

curveArray GraMath::To3d(curveArray cvs, geoMatrix3d mat)
{
	curveArray cvsVer;
	for each (auto & cv in cvs)
	{
		cvsVer.push_back(To3d(cv, mat));
	}

	return cvsVer;
}

CurveDataAlum GraMath::To3d(CurveDataAlum cv, geoMatrix3d mat)
{
	return CurveDataAlum(To3d(cv.ps, mat), To3d(cv.pe, mat));
}

QGEO::geoPoint3d GraMath::To3d(QGEO::geoPoint3d pt, geoMatrix3d mat)
{
	QGEO::geoPoint3d pt3d = pt;

	geoMatrix3d nimat;
	if (mat.inverse(nimat, 0.000001))
	{
		pt3d.TransformBy(nimat);

		//GraMath::ValidPoint(pt3d);
	}

	return pt3d;
}

void GraMath::LineSplitedByLoops(CurveDataAlum cv, std::vector<ARCHIOUT_LOOPDATA_ALUM> lps,
	std::vector<CurveDataAlum>& cvsOutLoop,
	std::vector<CurveDataAlum>& cvsInLoop,
	std::vector<CurveDataAlum>& cvsOnLoop)
{

	loopArray lpsValid;
	for each (auto & lp in lps)
	{
		//if (IsCurveInteredLoopByExtend(cv, lp))
		//{
		lpsValid.push_back(lp);
		//}
	}


	// 统计所有分割线条
	std::vector<CurveDataAlum> cvsLoops;
	for each (auto & lp in lpsValid)
	{
		for each (auto & cv in lp.m_curves)
		{
			cvsLoops.push_back(cv);
		}
	}


	std::vector<CurveDataAlum> cvsCutted = LineSplitedByCurves(cv, cvsLoops);
	for each (auto & cvCut in cvsCutted)
	{
		if (GraMath::IsPointInOneOfLoops(cvCut.pm, lps, 0))
		{
			cvsInLoop.push_back(cvCut);
		}
		else
		{

			bool bOnSide = false;
			for each (auto & cvLoop in cvsLoops)
			{
				if (GraMath::IsPointOnLine2d(cvCut.pm, cvLoop))
				{
					cvsOnLoop.push_back(cvCut);
					bOnSide = true;
					break;
				}
			}

			if (!bOnSide)
			{
				cvsOutLoop.push_back(cvCut);
			}
		}
	}
}

bool GraMath::IsPointOnLine2d(QGEO::geoPoint3d pt, CurveDataAlum cv, double dte)
{
	return IsPointOnLine2d(cv.ps, cv.pe, pt, dte);
}

double GraMath::disToLine2d(QGEO::geoPoint3d pt, QGEO::geoPoint3d p1, QGEO::geoPoint3d p2)
{
	geoLine l1 = geoLine(p1, p2);

	QGEO::geoPoint3d per;
	l1.GetPerpPt(pt, per);

	return GraMath::Distance2d(pt, per);
}


bool GraMath::IsPointOnLine2d(QGEO::geoPoint3d ps2, QGEO::geoPoint3d pe2, QGEO::geoPoint3d pt, double de)
{
	pe2.z = ps2.z;
	pt.z = ps2.z;

	double dis = GraMath::disToLine2d(pt, ps2, pe2);
	if (dis < 10)
	{
		if (GraMath::IsPointPerInLine(ps2, pe2, pt, de))
		{
			return true;
		}
	}

	return false;
}

bool GraMath::IsPointPerInLine(QGEO::geoPoint3d ps2, QGEO::geoPoint3d pe2, QGEO::geoPoint3d pt, double de)
{
	QGEO::geoPoint3d per = GraMath::PerToLine(pt, ps2, pe2);
	return IsPointInLine(per, ps2, pe2, de);
}


bool GraMath::IsPointInLine(QGEO::geoPoint3d& pt, QGEO::geoPoint3d ps2, QGEO::geoPoint3d pe2, double de)
{
	bool bOK = false;
	if (pt.DistanceTo(ps2) < de) bOK = true;
	if (pt.DistanceTo(pe2) < de) bOK = true;

	QGEO::geoVector3d v1 = (ps2 - pt).Normalize();
	QGEO::geoVector3d v2 = (pe2 - pt).Normalize();

	double ang0 = v1.AngleTo(v2);
	if (fabs(ang0) > .1) bOK = true;

	return bOK;
}

std::vector<QGEO::geoPoint3d> GraMath::OffsetPoints(std::vector<QGEO::geoPoint3d> pts, double wid, bool bInside)
{
	if (fabs(wid) < 1) return pts;

	std::vector<QGEO::geoPoint3d> ptsOff;

	if (pts.size() == 0) return ptsOff;

	std::vector<double> widList;


	int npt = pts.size();
	for (int j = 0; j < npt; j++)
	{
		widList.push_back(wid);
	}

	GraMath::OffsetPoly(widList, pts, ptsOff, bInside);

	return ptsOff;
}


void GraMath::OffsetPoly(std::vector<double> widList, std::vector<QGEO::geoPoint3d> pts, std::vector<QGEO::geoPoint3d>& ptsNew, bool bInside)
{
	GraMath::MakePtsClosed(pts);

	ARCHIOUT_LOOPDATA_ALUM lp(pts);



	int ncv = pts.size() - 1;

	int nwid = widList.size();
	if (nwid == 0) return;
	if (nwid < ncv)
	{
		for (int i = 0; i < ncv - nwid; i++)
		{
			widList.push_back(widList[0]);
		}
	}

	ARCHIOUT_LOOPDATA_ALUM wallFaceLoop = ARCHIOUT_LOOPDATA_ALUM();

	for (int i = 0; i < ncv; i++)
	{
		int iPre = (i - 1 + ncv) % ncv;
		int iNex = (i + 1) % ncv;
		int iNex2 = (i + 2) % ncv;

		CurveDataAlum cv = CurveDataAlum(pts[i], pts[iNex]);
		CurveDataAlum cvPre = CurveDataAlum(pts[iPre], pts[i]);
		CurveDataAlum cvNex = CurveDataAlum(pts[iNex], pts[iNex2]);

		QGEO::geoVector3d vecInner = GraMath::vec2dInner(cv, lp);
		QGEO::geoVector3d vecInnerPre = GraMath::vec2dInner(cvPre, lp);
		QGEO::geoVector3d vecInnerNex = GraMath::vec2dInner(cvNex, lp);

		QGEO::geoVector3d vecOff = (bInside) ? vecInner : -vecInner;
		QGEO::geoVector3d vecPreOff = (bInside) ? vecInnerPre : -vecInnerPre;
		QGEO::geoVector3d vecNexOff = (bInside) ? vecInnerNex : -vecInnerNex;

		CurveDataAlum cvOff = OffsetCurve(cv, vecOff, widList[i]);
		CurveDataAlum cvPreOff = OffsetCurve(cvPre, vecPreOff, widList[iPre]);
		CurveDataAlum cvNexOff = OffsetCurve(cvNex, vecNexOff, widList[iNex]);

		if (i == ncv - 1)
		{
			CurveTrimedByTwoSide(wallFaceLoop.m_curves, cvOff, cvPreOff, cvNexOff, true);
		}
		else
		{
			CurveTrimedByTwoSide(wallFaceLoop.m_curves, cvOff, cvPreOff, cvNexOff, false);
		}
	}

	GraMath::CurvesToPts(wallFaceLoop.m_curves, ptsNew);
}


void GraMath::CurvesToPts(std::vector<CurveDataAlum> cvs, std::vector<QGEO::geoPoint3d>& pts)
{
	pts.clear();

	int len = cvs.size();
	for (int i = 0; i < len; i++)
	{
		GraMath::ValidPoint(cvs[i].ps);

		pts.push_back(cvs[i].ps);
		if (i == len - 1)
		{
			GraMath::ValidPoint(cvs[i].pe);

			pts.push_back(cvs[i].pe);
		}
	}
}

void GraMath::ValidPoint(QGEO::geoPoint3d& pt)
{
	double xx = pt.x;
	double yy = pt.y;
	double zz = pt.z;

	ValidData(xx);
	ValidData(yy);
	ValidData(zz);

	pt.x = GetInt(xx);
	pt.y = GetInt(yy);
	pt.z = GetInt(zz);
}



void GraMath::ValidData(double& da)
{
	int dte = 10000;


	double dda = abs(da);

	if (dda > 1000000)
	{
		dte = 1;
	}
	else if (dda > 100000)
	{
		dte = 10;
	}
	else if (dda > 10000)
	{
		dte = 100;
	}
	else if (dda > 1000)
	{
		dte = 1000;
	}
	else if (dda > 100)
	{
		dte = 10000;
	}


	double xx1 = GetInt(da * dte) / (dte + 0.0);
	da = xx1;
}


QGEO::geoVector3d GraMath::vec2dInner(CurveDataAlum cvSide, ARCHIOUT_LOOPDATA_ALUM loop, double dteInner)
{
	QGEO::geoVector3d vecLeft = GraMath::vec2dLeft(cvSide);
	QGEO::geoPoint3d pTest = cvSide.pm + vecLeft * dteInner;

	if (GraMath::IsPointInPolyGon2d(pTest, loop))
	{
		return vecLeft;
	}
	else return -vecLeft;
}

BOOL GraMath::IsPointInPolyGon2d(QGEO::geoPoint3d pt, ARCHIOUT_LOOPDATA_ALUM lp)
{
	return IsPointInPolyGon2d(lp.m_pts, pt);
}


//************************************
// 函数名: PublicFunc::IsPointInPolyGon2d
// 返回:   BOOL
//
// 参数:   std::vector<QGEO::geoPoint3d> pts
// 参数:   QGEO::geoPoint3d pt
// 
// 功能:   判断点是否在环内
// 编写: 
//************************************
BOOL GraMath::IsPointInPolyGon2d(std::vector<QGEO::geoPoint3d> pts3, QGEO::geoPoint3d pt3d)
{
	int npts = pts3.size();
	if (npts < 3) return FALSE;
	UINT isInside =
		PointInsidePolygon(pts3, pt3d);

	if (isInside == 1) return TRUE;
	return FALSE;
}


UINT GraMath::PointInsidePolygon(std::vector<QGEO::geoPoint3d> points, const QGEO::geoPoint3d& pt, double EPSILON)
{
	int numPoints = points.size();
	// result = 0  点在多边形外
	// result = 1  点在多边形内
	// result = 2  点在多边形边上

	if (numPoints <= 3) return 0;

	double xMin = points[0].x;
	double yMin = points[0].y;
	double xMax = points[0].x;
	double yMax = points[0].y;
	UINT i;
	for (i = 1; i < numPoints; i++)
	{
		if (points[i].x < xMin) xMin = points[i].x;
		if (points[i].y < yMin) yMin = points[i].y;
		if (points[i].x > xMax) xMax = points[i].x;
		if (points[i].y > yMax) yMax = points[i].y;
	}

	if (pt.x < xMin - EPSILON
		|| pt.x > xMax + EPSILON
		|| pt.y < yMin - EPSILON
		|| pt.y > yMax + EPSILON)
	{
		return 0; // outside of polygon
	}

	double dd = 0.0;
	for (i = 0; i < numPoints; i++)
	{
		dd = (fabs(pt.x - points[i].x) + fabs(pt.y - points[i].y));
		if (dd < (2.0 * EPSILON))
		{
			return 2;
		}
	}

	double DY, DDY, DDP, DII, DJJ, XPP, DX, DL, DKK;
	int num;
	int flag;
	int II, L;

	num = 0;
	flag = 0;
	UINT j = 0;
	for (i = 0; i < numPoints; i++)
	{
		if (i == numPoints - 1)
			j = 0;
		else
			j = i + 1;
		if (
			(points[i].y > pt.y + EPSILON && points[j].y > pt.y + EPSILON)
			|| (points[i].y < pt.y - EPSILON && points[j].y < pt.y - EPSILON)
			|| (points[i].x < pt.x - EPSILON && points[j].x < pt.x - EPSILON)) continue;

		DY = points[j].y - points[i].y;
		DDY = fabs(DY);
		DDP = fabs(pt.y - points[i].y);
		if (DDY < EPSILON && DDP < EPSILON)
		{
			xMin = min(points[i].x, points[j].x) - EPSILON;
			xMax = max(points[i].x, points[j].x) + EPSILON;
			if (pt.x >= xMin && pt.x < xMax)
			{
				return 2;
			}
			continue;
		}
		else if (DDY < EPSILON)
		{
			continue;
		}
		DX = points[j].x - points[i].x;
		DL = -DX * points[i].y + DY * points[i].x;
		XPP = (DX * pt.y + DL) / DY;
		dd = XPP - pt.x;

		if (fabs(dd) < EPSILON)
		{
			return 2;
		}
		else if (XPP > pt.x)
		{
			if (flag > 0)
			{
				DII = DJJ;
			}
			else
			{
				DII = sqrt((XPP - points[i].x) * (XPP - points[i].x) + (pt.y - points[i].y) * (pt.y - points[i].y));
				flag = 1;
			}
			DJJ = sqrt((XPP - points[j].x) * (XPP - points[j].x) + (pt.y - points[j].y) * (pt.y - points[j].y));
			if (DJJ < EPSILON)
			{
				num = num + 1;
			}
			else if (DII < EPSILON)
			{
				II = i;
			MARK20:
				if (II == 0)
					L = numPoints - 1;
				else
					L = II - 1;
				if (fabs(points[L].y - points[i].y) < EPSILON)
				{
					II = L;
					goto MARK20;
				}
				else if (
					(points[L].y - EPSILON < pt.y
						&& points[j].y - EPSILON < pt.y)
					|| (points[L].y + EPSILON > pt.y
						&& points[j].y + EPSILON > pt.y))
				{
					num = num - 1;
				}
				continue;
			}
			else
			{
				DKK = (XPP - points[j].x) * (points[i].x - XPP)
					+ (pt.y - points[j].y) * (points[i].y - pt.y);
				if (DKK < 0.0) continue;
				num = num + 1;
			}
		}
	}

	return num % 2;
}


bool GraMath::IsPointInOneOfLoops(QGEO::geoPoint3d pm, std::vector<ARCHIOUT_LOOPDATA_ALUM> bounds, int offsetOutdis)
{

	if (bounds.empty()) return false;

	if (GraMath::IsFaceHorizon(bounds[0]))
	{
		int len = bounds.size();
		for (int i = 0; i < len; i++)
		{
			std::vector<QGEO::geoPoint3d> ptsOff = GraMath::OffsetPoints(bounds[i].m_pts, offsetOutdis, false);

			if (GraMath::IsPointInPolyGon2d(ptsOff, pm))
			{
				return true;
			}
		}
	}
	else
	{
		geoMatrix3d mat = matrixFrom(bounds[0]);

		QGEO::geoPoint3d pm2d = GraMath::To2d(pm, mat);
		loopArray bounds2d = GraMath::To2d(bounds, mat);

		int len = bounds2d.size();
		for (int i = 0; i < len; i++)
		{
			std::vector<QGEO::geoPoint3d> ptsOff = GraMath::OffsetPoints(bounds2d[i].m_pts, offsetOutdis, false);

			if (GraMath::IsPointInPolyGon2d(ptsOff, pm2d))
			{
				return true;
			}
		}
	}

	return false;
}

geoMatrix3d GraMath::matrixFrom(ARCHIOUT_LOOPDATA_ALUM lp)
{
	return matrixFrom(lp.m_pts);
}

bool GraMath::IsFaceHorizon(ARCHIOUT_LOOPDATA_ALUM lp)
{
	return IsFaceHorizon(lp.m_pts);
}

bool GraMath::IsFaceHorizon(const std::vector<QGEO::geoPoint3d>& vertexs)
{
	QGEO::geoVector3d vecNor = normalPolyGon(vertexs);
	if (fabs(fabs(vecNor.z) - 1) < 0.1) return true;
	return false;
}


void GraMath::SetCurvesZ(curveArray& cvs, double z)
{
	int ncvs = cvs.size();
	for (int i = 0; i < ncvs; i++)
	{
		SetCurveZ(cvs[i], z);
	}
}


BOOL GraMath::IsTwoLineParrallel(QGEO::geoPoint3d p1, QGEO::geoPoint3d p2, QGEO::geoPoint3d t1, QGEO::geoPoint3d t2, double dte)
{
	QGEO::geoVector3d v1 = (p1 - p2).Normalize();
	QGEO::geoVector3d v2 = (t1 - t2).Normalize();

	geoVector2d v1d(v1.x, v1.y);
	geoVector2d v2d(v2.x, v2.y);

	double diff = fabs(v1d.x * v2d.y - v1d.y * v2d.x);

	if (diff < dte) return TRUE;

	return FALSE;
}

bool GraMath::IsTwoLineParrallel(CurveDataAlum cv1, CurveDataAlum cv2, double dte)
{
	return IsTwoLineParrallel(cv1.ps, cv1.pe, cv2.ps, cv2.pe, dte);
}

std::vector<CurveDataAlum> GraMath::LineSplitedByCurves(CurveDataAlum cv, std::vector<CurveDataAlum> cvsCut)
{
	std::vector<CurveDataAlum> cvsAll, curvesResult, cvsSplited;


	GraMath::SetCurvesZ(cvsCut, cv.ps.z);

	CurveDataAlum cvm(cv.ps, cv.pe, DWORD_PTR(-1001));
	cvsAll.push_back(cvm);

	//std::vector<CurveDataAlum> cvsCutValid;
	for (auto& cvCut : cvsCut)
	{
		if (GraMath::IsTwoLineParrallel(cv, cvCut, 0.01)) continue;
		cvsAll.push_back(cvCut);
	}


	//GraMath::vecAppend(cvsAll, cvsCut);
	SplitCurves(cvsAll, curvesResult, 1.0);

	for each (auto cvSplit in curvesResult)
	{
		SetCurveZ(cvSplit, cv.ps.z);
		int udata = (int)cvSplit.userData;
		if (udata == -1001)
		{
			cvsSplited.push_back(cvSplit);
		}
	}

	return cvsSplited;
}

ARCHIOUT_LOOPDATA_ALUM GraMath::GetRectLoopVert(const QGEO::geoPoint3d& ps, const QGEO::geoPoint3d& pe, double zmin, double zmax)
{
	double xs = ps.x;
	double ys = ps.y;
	double xe = pe.x;
	double ye = pe.y;

	QGEO::geoPoint3dArray pts =
	{
		QGEO::geoPoint3d(xs, ys, zmin),
		QGEO::geoPoint3d(xe, ye, zmin),
		QGEO::geoPoint3d(xe, ye, zmax),
		QGEO::geoPoint3d(xs, ys, zmax),
		QGEO::geoPoint3d(xs, ys, zmin),
	};

	return ARCHIOUT_LOOPDATA_ALUM(pts);
}

ARCHIOUT_LOOPDATA_ALUM GraMath::GetRectLoopVert(CurveDataAlum cv, double zmin, double zmax)
{
	return GetRectLoopVert(cv.ps, cv.pe, zmin, zmax);
}


void GraMath::SplitCurves(std::vector<CurveDataAlum> curves, std::vector<CurveDataAlum>& curvesResult, double tol_point_eq)
{
	if (curves.size() == 0) return;

	//curvesResult = curves;

	ISearchLoop searchLoop;
	QGEO::geoPoint3d pointA, pointB;
	for (auto curve : curves)
	{
		AddCurveToFindLoop(searchLoop, curve);
	}

	ISearchLoopResult* result = NULL;
	searchLoop.GetSplitCurveResult(result, tol_point_eq);

	if (result == NULL)
	{
		curvesResult = curves;
		return;
	}

	tvector<GeaCurve2d*> arrCV;
	tvector<DWORD_PTR> userDatas;

	result->getUserDataCurve(arrCV, userDatas);
	result->getCurve(arrCV);
	int nline = arrCV.size();
	int ndata = userDatas.size();

	for (int i = 0; i < nline; i++)
	{
		if (arrCV[i] == nullptr) continue;

		Pnt2d ps = arrCV[i]->getEndPoint_(0);
		Pnt2d pe = arrCV[i]->getEndPoint_(1);

		DWORD_PTR ud = (i < ndata) ? userDatas[i] : 0;

		QGEO::geoPoint3d pps(ps.x, ps.y, 0);
		QGEO::geoPoint3d ppe(pe.x, pe.y, 0);
		if (arrCV[i]->classType_() == GeaCurve2d::CV_LINESEG)
		{
			CurveDataAlum cv(pps, ppe, ud);
			curvesResult.push_back(cv);
		}
		else if (arrCV[i]->classType_() == GeaCurve2d::CV_ARC)
		{
			ARC2D_DATA data;
			arrCV[i]->getArcData(data);
			QGEO::geoPoint3d pc(data.m_center.x, data.m_center.y, 0);
			bool bAntiWise = (data.m_sweep > M_PI) ? false : true;

			CurveDataAlum cvArc;
			cvArc.isArc = true;
			cvArc.userData = ud;

			cvArc.ps = pps;
			cvArc.pe = ppe;
			curvesResult.push_back(cvArc);
		}
	}

	result->release();
}

void GraMath::AddCurveToFindLoop(ISearchLoop& searchLoop, CurveDataAlum cv)
{
	if (cv.isArc)
	{
		geoArc arc;
		arc.set3PointArc(cv.ps, cv.pm, cv.pe);
		AddArcToFindLoop(searchLoop, arc, cv.userData);
	}
	else
	{
		searchLoop.addLine(Pnt2d(cv.ps.x, cv.ps.y), Pnt2d(cv.pe.x, cv.pe.y), cv.userData);
	}
}



void GraMath::AddArcToFindLoop(ISearchLoop& searchLoop, geoArc arc, DWORD_PTR userData)
{
	QGEO::geoPoint3d pc = arc.center();
	double rad = arc.radius();
	QGEO::geoPoint3d ps = arc.startPoint();
	QGEO::geoPoint3d pe = arc.endPoint();
	QGEO::geoVector3d vs = (ps - pc).Normalize();
	QGEO::geoVector3d ve = (pe - pc).Normalize();
	QGEO::geoVector3d vNor = arc.normal();
	QGEO::geoVector3d vecStartTan = vs;

	QGEO::geoVector3d gvs(vs.x, vs.y, 0);
	QGEO::geoVector3d gve(ve.x, ve.y, 0);

	double as = AngleOfDVec3d(gvs) / ARADK;
	double ae = AngleOfDVec3d(gve) / ARADK;

	if (fabs(vNor.z) > 0)
	{
		if (as > ae) ae += 360;
	}

	int nSplit = (int)((ae - as) / (15 + 0.0));

	std::vector<QGEO::geoPoint3d> ptsArc;
	ArcToPtList(pc, rad, as, ae, ptsArc, nSplit);

	int len = ptsArc.size();
	for (int i = 0; i < len - 1; i++)
	{
		QGEO::geoPoint3d s1 = ptsArc[i];
		QGEO::geoPoint3d s2 = ptsArc[i + 1];

		searchLoop.addLine(Pnt2d(s1.x, s1.y), Pnt2d(s2.x, s2.y), userData);

		if (i == len - 2)
		{
			searchLoop.addLine(Pnt2d(s2.x, s2.y), Pnt2d(pe.x, pe.y), userData);
		}
	}

	return;

	double sweepAngle2 = gvs.AngleTo(gve);
	double sweepAngle = vs.AngleTo(ve);
	if (fabs(vNor.z) < 0)
	{
		for (double ang = as; ang < ae; ang += 10)
		{
		}

		sweepAngle = 2 * M_PI - sweepAngle;
		vecStartTan = ve;
	}
	else
	{
		if (sweepAngle < 0)
		{
			sweepAngle = 2 * M_PI + sweepAngle;
		}
	}

	sweepAngle = 1.2 * M_PI /*- sweepAngle*/;

	vecStartTan = ve;

	searchLoop.addArc(Pnt2d(pc.x, pc.y), rad, Vec2d(vecStartTan.x, vecStartTan.y), sweepAngle, userData);
}


void GraMath::ArcToPtList(QGEO::geoPoint3d pc, float r, float as, float ae, std::vector<QGEO::geoPoint3d>& pts, int nSplit)
{
	pts.clear();
	if (as > ae)
	{
		ae += 360;
	}

	double angadd = (ae - as) / (nSplit + 0.0);
	for (size_t i = 0; i < nSplit + 1; i++)
	{
		double ang = as + i * angadd;
		QGEO::geoVector3d vec(cos(ang * ARADK) * r, sin(ang * ARADK) * r, 0);
		QGEO::geoPoint3d pp = FromSumOf(pc, vec);
		pts.push_back(pp);
	}
}

CurveDataAlum GraMath::To2d(CurveDataAlum cv, geoMatrix3d mat)
{
	return CurveDataAlum(To2d(cv.ps, mat), To2d(cv.pe, mat));
}

loopArray GraMath::To2d(loopArray lps, geoMatrix3d mat)
{
	loopArray lps2;
	for each (auto lp in lps)
	{
		lps2.push_back(To2d(lp, mat));
	}

	return lps2;
}

ARCHIOUT_LOOPDATA_ALUM GraMath::To2d(ARCHIOUT_LOOPDATA_ALUM lp, geoMatrix3d mat)
{
	QGEO::geoPoint3dArray ptsNew;
	for each (auto pp in lp.m_pts)
	{
		ptsNew.push_back(To2d(pp, mat));
	}
	return ARCHIOUT_LOOPDATA_ALUM(ptsNew);
}

QGEO::geoPoint3d GraMath::To2d(QGEO::geoPoint3d pt, geoMatrix3d mat)
{
	QGEO::geoPoint3d pt2d = pt;
	pt2d.TransformBy(mat);

	//GraMath::ValidPoint(pt2d);

	return pt2d;
}

geoMatrix3d GraMath::matrixFrom(std::vector<QGEO::geoPoint3d> points)
{
	if (points.size() < 3) return geoMatrix3d();

	geoMatrix3d mat;
	mat.setToWorldToPlane(normalPolyGon(points), points[0]);

	return mat;
}

BOOL GraMath::IsFileExist(CString strFileName)
{
	DWORD dwAttribute = GetFileAttributes(strFileName);
	if (dwAttribute == INVALID_FILE_ATTRIBUTES || dwAttribute & FILE_ATTRIBUTE_DIRECTORY)
		return FALSE;
	return TRUE;
}


CString GraMath::GetFileName(const CString& filePathName)
{
	CString firstFileName = _T("");
	firstFileName = filePathName.Right(filePathName.GetLength() - filePathName.ReverseFind('\\') - 1);
	return firstFileName;
}

QGEO::geoPoint3d GraMath::GetCentroid(const std::vector<QGEO::geoPoint3d>& vertexs)
{
	if (IsFaceHorizon(vertexs))
	{
		return GraMath::GetCentroid2d(vertexs);
	}
	else
	{
		geoMatrix3d mat = GraMath::matrixFrom(vertexs);

		ARCHIOUT_LOOPDATA_ALUM lp(vertexs), lpPlan;
		lpPlan = GraMath::To2d(lp, mat);
		QGEO::geoPoint3d pcPlan = GraMath::GetCentroid2d(lpPlan.m_pts);
		return To3d(pcPlan, mat);
	}
}


QGEO::geoPoint3d GraMath::GetCentroid2d(const std::vector<QGEO::geoPoint3d>& vertexs)
{
	QGEO::geoPoint3d centroid(0, 0, 0);

	int vertexCount = vertexs.size();

	if (vertexCount == 0) return centroid;

	int z0 = vertexs[0].z;


	double signedArea = 0.0;
	double x0 = 0.0; // 当前顶点 X
	double y0 = 0.0; // 当前顶点 Y
	double x1 = 0.0; // 下一顶点 X
	double y1 = 0.0; // 下一顶点 Y
	double a = 0.0;  // Partial signed area

	// 遍历所有顶点（首尾重合，不包含结尾）
	int i = 0;
	for (i = 0; i < vertexCount - 1; ++i)
	{
		x0 = vertexs[i].x;
		y0 = vertexs[i].y;
		x1 = vertexs[i + 1].x;
		y1 = vertexs[i + 1].y;
		a = x0 * y1 - x1 * y0;
		signedArea += a;
		centroid.x += (x0 + x1) * a;
		centroid.y += (y0 + y1) * a;
	}

	//单独处理最后一个顶点，以避免在每次迭代中执行复杂的运算
	x0 = vertexs[i].x;
	y0 = vertexs[i].y;
	x1 = vertexs[0].x;
	y1 = vertexs[0].y;
	a = x0 * y1 - x1 * y0;
	signedArea += a;
	centroid.x += (x0 + x1) * a;
	centroid.y += (y0 + y1) * a;

	signedArea *= 0.5;
	centroid.x /= (6.0 * signedArea);
	centroid.y /= (6.0 * signedArea);
	centroid.z = z0;
	return centroid;
}

BOOL GraMath::IsStrHasStr(CString cs, CString csKey)
{
	int id1 = 0, id2 = 0;
	return IsStrHasStr(cs, csKey, id1, id2);
}


BOOL GraMath::IsStrHasStr(CString cs, CString cs1, int& idhead, int& idtail)
{
	idhead = 0;
	idtail = 0;

	int istart = cs.Find(cs1);
	int len = cs1.GetLength();

	if (istart != -1)
	{
		idhead = istart;
		idtail = idhead + len;
		return TRUE;
	}

	return FALSE;
}

void GraMath::SetDrawInDrag(bool bDrag)
{
	g_bInDrag = bDrag;
}

bool GraMath::GetDrawInDrag()
{
	return g_bInDrag;
}


