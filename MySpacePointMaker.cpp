// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "MySpacePointMaker.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
//#include "ActsExamples/Utilities/GroupBy.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <list>
#include "TMath.h"
//#include "tracker_config.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"

ActsExamples::MySpacePointMaker::MySpacePointMaker(Config cfg, Acts::Logging::Level lvl)
  : IAlgorithm("MySpacePointMaker", lvl), m_cfg(std::move(cfg)) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);

  auto spConstructor =
      [](const Acts::Vector3& pos, std::optional<double> t,
         const Acts::Vector2& cov, std::optional<double> varT,
         boost::container::static_vector<Acts::SourceLink, 2> slinks)
      -> SimSpacePoint {
    return SimSpacePoint(pos, t, cov[0], cov[1], varT, std::move(slinks));
  };

  auto spBuilderConfig = Acts::SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = m_cfg.trackingGeometry;
  m_slSurfaceAccessor.emplace(IndexSourceLink::SurfaceAccessor{*m_cfg.trackingGeometry});
  spBuilderConfig.slSurfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&m_slSurfaceAccessor.value());
  m_spacePointBuilder = Acts::SpacePointBuilder<SimSpacePoint>(spBuilderConfig, spConstructor, Acts::getDefaultLogger("SpacePointBuilder", lvl));
}

ActsExamples::ProcessCode ActsExamples::MySpacePointMaker::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  // function to access measurement parameters using source links
  auto accessor = [&measurements](Acts::SourceLink slink) {
    const auto islink = slink.get<IndexSourceLink>();
    const ConstVariableBoundMeasurementProxy meas = measurements.getMeasurement(islink.index());
    return std::make_pair(meas.fullParameters(), meas.fullCovariance());
  };

  // fill front and back strip vectors of pairs(source link, pair<stripEnd1, stripEnd2>)
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> frontStrips;
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> backStrips;
  std::vector<Acts::SourceLink> twoDimMeasurements;

  int shift = 3; // FIXME
  int nStations = 5; // TODO read from config
  std::list<std::vector<ActsExamples::IndexSourceLink>> candidates[nStations];

  double rMax = 1300; // FIXME Extract from surface dimensions
  for (auto& isl : measurements.orderedIndices()) {
    const auto geoId = isl.geometryId();
    const auto volumeId = geoId.volume();
    const auto layerId = geoId.layer();
    int iStation = (layerId-shift)/7;
    int iLayer = (layerId-shift)%7;

    ACTS_DEBUG("volumeId= " << volumeId << " layerId=" << layerId);
    Acts::SourceLink slink{isl};

    if ((layerId-shift)%7==3 || layerId==2 || layerId==38) {
    // if (layerType[layerId - shift]==2) {
      twoDimMeasurements.emplace_back(slink);
    } else if (0){
      const auto [par, cov] = accessor(slink);
      const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink);
      // TODO: more realistic strip dimensions including inner radii and half-station splitting
      auto gpos1 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0],-rMax), Acts::Vector3());
      auto gpos2 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0], rMax), Acts::Vector3());
      if (layerId%4==2) {
        frontStrips.emplace_back(std::make_pair(slink, std::make_pair(gpos1, gpos2)));
      } else if (layerId%4==0) {
        backStrips.emplace_back(std::make_pair(slink, std::make_pair(gpos1, gpos2)));
      }
    } else {
      int iStraw = geoId.sensitive();
      // construct sp candidates as vectors of index source links from the same station corresponding to
      // straws with similar straw ids (~similar angles)
      for (auto& candidate : candidates[iStation]){
        bool isCompatibleCandidate = 1;
        for (auto& isl2 : candidate){
          const auto geoId2 = isl2.geometryId();
          int iLayer2 = (geoId2.layer()-shift)%7;
          int iStraw2 = geoId2.sensitive();
          int layerDif = iLayer - iLayer2;
          int strawDif = iStraw - iStraw2;
          if (layerDif<=0) { isCompatibleCandidate = 0; break ; } // should be imposible by construction
          if (layerDif==4) { // same direction
            if (strawDif<-1 || strawDif>0) { isCompatibleCandidate = 0; break ; }
          } else if (layerDif==1){
            if ((iLayer2==0 || iLayer2==4) && (strawDif<-2 || strawDif>3)) { isCompatibleCandidate = 0; break ; };
            if ((iLayer2==1 || iLayer2==5) && (strawDif<-6 || strawDif>4)) { isCompatibleCandidate = 0; break ; };
          } else if (layerDif==2){
            if ((iLayer2==0 || iLayer2==4) && (strawDif<-4 || strawDif>2)) { isCompatibleCandidate = 0; break ; };
          }
        }
        if (!isCompatibleCandidate) continue;
        candidate.push_back(isl);
        ACTS_DEBUG("Adding new measurement to existing candidate");
      }
      if (iLayer==0 || iLayer==1 || iLayer==2 || iLayer==4){
        std::vector<ActsExamples::IndexSourceLink> new_candidate;
        new_candidate.push_back(isl);
        candidates[iStation].push_back(new_candidate);
        ACTS_DEBUG("Adding new candidate");
      }
    }
  }

  for (int iStation=0;iStation<nStations;iStation++) {
    ACTS_DEBUG("Station " << iStation << ": number of candidates " << candidates[iStation].size());    
    // for (auto candidate : candidates[iStation]) {
    for (auto it = candidates[iStation].begin(); it != candidates[iStation].end();) {
      auto& candidate = *it;
      ACTS_DEBUG("  candidate.size=" << candidate.size());
      for (auto& isl : candidate) {
        ACTS_DEBUG("   layer=" << isl.geometryId().layer() << " straw=" << isl.geometryId().sensitive());
      }
      if (candidate.size()<3) {
        ACTS_DEBUG(  "erasing...");
        it = candidates[iStation].erase(it);
      } else {
        it++;
      }
    }
  }
  ACTS_DEBUG("making space points from straws");

  std::vector<double> c;
  std::vector<double> s;
  std::vector<double> z;
  std::vector<double> g;      
  std::vector<double> d;

  for (int iStation=0;iStation<nStations;iStation++) {
    ACTS_DEBUG("Station " << iStation << ": number of filtered candidates " << candidates[iStation].size());    
    for (auto& candidate : candidates[iStation]) {
      ACTS_VERBOSE("  candidate.size=" << candidate.size());
      int n = candidate.size();
      c.resize(n);
      s.resize(n);
      z.resize(n);
      g.resize(n);
      d.resize(n);
      for (size_t i=0; i<n; i++){
        auto& isl = candidate[i];
        const auto geoId = isl.geometryId();
        const auto layerId = geoId.layer();
        const auto strawId = geoId.sensitive();
        Acts::SourceLink slink{isl};
        const auto [par, cov] = accessor(slink);
        // tube center shifted by the measured distance to wire
        const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink);
        auto xyz = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0],0), Acts::Vector3());
        auto rot = surface->transform(ctx.geoContext).rotation();
        double x = xyz[0];
        double y = xyz[1];
        double cosp = rot(0,1);
        double sinp = rot(0,0);
        z[i] = xyz[2];
        s[i] = z[i]*sinp;
        c[i] = z[i]*cosp;
        g[i] = -x*sinp + y*cosp;
        d[i] = sqrt(cov(0,0));
      }
      double t,k,dt,dk;
      double chi2 = analytic(s, c, g, d, t, k, dt, dk);
      ACTS_DEBUG("  n=" << n << " t=" << t <<" k=" << k << " chi2/ndf=" << chi2/(n-2));
      double tx,ty,p;
      double chi2helix = helix(z, s, c, g, d, tx, ty, p);
//      ACTS_DEBUG("  n=" << n << " tx=" << tx <<" ty=" << ty << " chi2/ndf=" << chi2helix/(n-2));
      ACTS_DEBUG("  n=" << n << " t=" << sqrt(tx*tx+ty*ty) <<" k=" << ty/tx << " chi2/ndf=" << (n>3 ? chi2helix/(n-3) : chi2helix));
    }
  }



  ACTS_VERBOSE("making strip pairs:" << "back: " << frontStrips.size() << " front: " << backStrips.size());

  // make space points from strips
  SimSpacePointContainer spacePoints;
  for (auto& fstrip : frontStrips) {
    float fz = fstrip.second.second[2];
    for (auto& bstrip : backStrips) {
      float bz = bstrip.second.second[2];
      if (fabs(fz - bz) > 20) continue;
      std::vector<Acts::SourceLink> slinks = {fstrip.first, bstrip.first};
      auto strippair = std::make_pair(fstrip.second, bstrip.second);
      Acts::SpacePointBuilderOptions spOptStrips{strippair, accessor};
      m_spacePointBuilder.buildSpacePoint(ctx.geoContext, slinks, spOptStrips, std::back_inserter(spacePoints));
    }
  }

  // make space points from 2D measurements
  Acts::SpacePointBuilderOptions spOpt;
  spOpt.paramCovAccessor = accessor;

  for (auto& slink : twoDimMeasurements) {
    m_spacePointBuilder.buildSpacePoint(ctx.geoContext, {slink}, spOpt, std::back_inserter(spacePoints));
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}


double ActsExamples::MySpacePointMaker::analytic(std::vector<double> &a, std::vector<double> &cz, std::vector<double> &g, std::vector<double> &s, 
                                                 double &t, double &k, double &dt, double &dk, bool debug) const{
  
  int n = a.size();
  std::vector<double> b(n,0);
  for (int i=0;i<n;i++) b[i] = -cz[i];
  std::vector<double> aa(n,0);
  std::vector<double> bb(n,0);
  double A = 0;
  double B = 0;
  double sigA2 = 0;
  double sigB2 = 0;
  double covAB = 0;

  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      double ab = (a[i]*b[j]-b[i]*a[j]);
      aa[i]+=b[j]*ab;
      bb[i]+=a[j]*ab;
    }
    A+= g[i]*aa[i];
    B+= g[i]*bb[i];
    sigA2+=s[i]*s[i]*aa[i]*aa[i];
    sigB2+=s[i]*s[i]*bb[i]*bb[i];
    covAB+=s[i]*s[i]*aa[i]*bb[i];
  }

  k = - B/A;
  dk = sqrt(B*B*sigA2+A*A*sigB2-2*A*B*covAB)/A/A;

  std::vector<double> abk(n,0);
  std::vector<double> vt(n,0);
  std::vector<double> dtdk(n,0);
  std::vector<double> dkdc(n,0);
  std::vector<double> dtdc(n,0);
  double num = 0;
  double den = 0;
  double p1 = 0;
  double p2 = 0;
  for (int i=0;i<n;i++){
    abk[i] = a[i]+b[i]*k;
    num += g[i]*abk[i];
    den += abk[i]*abk[i];
    p1 += a[i]*b[i];
    p2 += b[i]*b[i];
  }
  double cosa = 1./sqrt(1+k*k);
  double sina = k*cosa;

  t = -1./cosa*num/den;

  for (int i=0; i<n; i++){
    vt[i]  = -1./cosa*abk[i]/den;
    dtdk[i] = vt[i]*(k/(1+k*k) + b[i]/abk[i]-2*(p1+p2*k)/den);
    dkdc[i] = (aa[i]*B - bb[i]*A)/A/A;
  }

  double dt2 = 0;
  for (int i=0; i<n; i++){
    dtdc[i] = vt[i];
    for (int j=0; j<n; j++){
      dtdc[i] += g[j]*dtdk[j]*dkdc[i];
    }
    dt2+=s[i]*s[i]*dtdc[i]*dtdc[i];
  }

  dt = sqrt(dt2);
  
  double chi2 = 0;
  for (int i=0;i<n;i++){
    double d = a[i]*t*cosa + b[i]*t*sina + g[i];
    chi2+=d*d/s[i]/s[i];
  }
  return chi2;
}

double ActsExamples::MySpacePointMaker::helix(
  std::vector<double> &z, std::vector<double> &s, std::vector<double> &c, std::vector<double> &g, 
  std::vector<double> &d, double &tx, double &ty, double &p, bool debug) const 
{
  int n = s.size();
  double sk[n]={0};
  double ck[n]={0};  
  auto fk = [&n,&s,&c,&z,&g,&tx,&ty,&sk,&ck](double k) {
    double ss = 0;
    double sc = 0;
    double cc = 0;
    double gs = 0;
    double gc = 0;
    for (int i=0;i<n;i++){
      sk[i] = (s[i] + k*z[i]*c[i]);
      ck[i] = (c[i] - k*z[i]*s[i]);
      ss += sk[i]*sk[i];
      sc += sk[i]*ck[i];
      cc += ck[i]*ck[i];
      gs += g[i]*sk[i];
      gc += g[i]*ck[i];
    }
    tx = (gc*sc - gs*cc)/(ss*cc - sc*sc);
    ty = (gc*ss - gs*sc)/(ss*cc - sc*sc);

    double sum = 0;
    for (int i=0;i<n;i++){
      sum+=(tx*sk[i]-ty*ck[i]+g[i])*(ty*z[i]*s[i]+tx*z[i]*c[i]);
    }
    return sum;
  };
  
  double f0 = fk(0);
  double dk = 1e-7;
  double dfdk = (fk(dk)-f0)/dk;
  double kk = -2*f0/dfdk;
  double kmin = kk>0 ? 0  : kk;
  double kmax = kk>0 ? kk :  0;
  // printf("%f %f %f %f\n",kmin, kmax, fk(kmin), fk(kmax));
  
  ROOT::Math::Functor1D functor(fk);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(functor, kmin, kmax);
  rf.Solve();
  // printf("k=%e iterations=%d\n",rf.Root(),rf.Iterations());

  p = 2*rf.Root();
  double chi2 = 0;
  for (int i=0;i<n;i++){
    double chi = (sk[i]*tx - ck[i]*ty + g[i])/d[i];
    chi2 += chi*chi;
  }
  return chi2;
}
