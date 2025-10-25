#ifndef tracker
#define tracker
#include "tracker_config.h"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/NavigationLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Utilities/Helpers.hpp"
//#include "Acts/Plugins/Geant4/Geant4Converters.hpp"

//#include "Geant4/G4Material.hh"
//#include "Geant4/G4NistManager.hh"


using Acts::UnitConstants::cm;
using namespace Acts::UnitLiterals;

class MyDetectorElement : public Acts::DetectorElementBase {
public:
  MyDetectorElement(std::shared_ptr<const Acts::Transform3> transform, std::shared_ptr<Acts::Surface> surface, double thickness)
      : Acts::DetectorElementBase(), mTransform(transform), mSurface(surface), mThickness(thickness)
  {
  }
  const Acts::Transform3 &transform(const Acts::GeometryContext &gctx) const { return *mTransform; }
  const Acts::Surface &surface() const { return *mSurface; }
  Acts::Surface &surface() { return *mSurface; }
  double thickness() const { return mThickness; }

private:
  std::shared_ptr<const Acts::Transform3> mTransform = nullptr;
  std::shared_ptr<Acts::Surface> mSurface = nullptr;
  double mThickness = 0.;
};

std::vector<std::shared_ptr<MyDetectorElement>> detectorStore;
Acts::GeometryContext gctx;

Acts::TrackingGeometry* CreateTrackingGeometry(bool addROC = 0, bool addFlange = 0){
  auto gctx = Acts::GeometryContext();

  auto silicon = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  auto vacuum = Acts::Material::fromMassDensity(1e+10_cm, 1e+10_cm, 1.0000, 1, 1e-10_g / 1_cm3);
  auto volumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(vacuum);

  std::vector<float> layerBinEdges;
  std::vector<std::pair<std::shared_ptr<const Acts::Layer>, Acts::Vector3>> layerOrderVec;
  // add first navigation layer
  Acts::Transform3 tr1 = Acts::Transform3::Identity();
  tr1.translate(Acts::Vector3(0, 0, positions[0]*cm - 20));
  const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMaxStation*cm, rMaxStation*cm);
  auto layer = Acts::PlaneLayer::create(tr1, pBounds, nullptr, 1_mm, nullptr, Acts::navigation);
  layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
  layerBinEdges.push_back(positions[0]*cm - 30);
  layerBinEdges.push_back(positions[0]*cm - 10);

  for (int i = 0; i < positions.size(); ++i) {
    double rmin = layerRMin[i] * cm;
    double rmax = layerRMax[i] * cm;
    layerBinEdges.push_back(i < positions.size()-1 ? (positions[i] + positions[i+1])/2.*cm : positions.back()*cm+ 1);

    if (layerType[i]==2) {
      const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rmax, rmax);
      Acts::Transform3 trafo = Acts::Transform3::Identity();
      trafo.translate(Acts::Vector3(0, 0, positions[i]*cm));
      auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds);
      double surfThickness = 0.0001 * cm;
      double layerThickness = 0.0001 * cm;
      auto matProp = Acts::MaterialSlab(silicon, surfThickness);
      auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);
      surface->assignSurfaceMaterial(std::move(surfaceMaterial));
      auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
      auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), layerThickness, nullptr, Acts::active);
      surface->associateLayer(*layer.get());
      layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
      auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, surfThickness);
      surface->assignDetectorElement(*detElement.get());
      detectorStore.push_back(std::move(detElement));
      continue;
    }

    double angleRot = 0;
    if (layerType[i]==5) angleRot = 7.*M_PI/180.;
    if (layerType[i]==6) angleRot =-7.*M_PI/180.;
    // create surface material
    double surfThickness = thickness * cm;
    double layerThickness = thickness * cm;
    auto matProp = Acts::MaterialSlab(silicon, surfThickness);
    auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);
    // create surface array
    // double rmin = 0.1 * cm;
    // double rmax = 150 * cm;
    double rc = 0.5 * (rmax + rmin);   
    double hl = 0.5 * (rmax - rmin);
    double hw = 0.2 * cm;
    const auto sBounds = std::make_shared<const Acts::RectangleBounds>(hw, hl);
    int nTubes = numberOfTubes[i];
    std::vector<std::shared_ptr<const Acts::Surface>> vSurfaces;
    for (int iTube = 0; iTube < nTubes; iTube++) {
      double phi = layerAngle[i] + 2*M_PI/nTubes*(iTube+0.5);
      double cosp = cos(phi);
      double sinp = sin(phi);
      auto stransform = Acts::Transform3::Identity();
      stransform.translate(Acts::Vector3(rc*cosp, rc*sinp, positions[i] * cm));
      stransform.rotate(Eigen::AngleAxisd(M_PI/2 + phi + angleRot, Acts::Vector3(0, 0, 1)));
      auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(stransform, sBounds);
      surface->assignSurfaceMaterial(surfaceMaterial);
      // create detector element
      auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(stransform), surface, surfThickness);
      surface->assignDetectorElement(*detElement);
      detectorStore.push_back(std::move(detElement));
      vSurfaces.push_back(std::move(surface));
    }
    auto trafo = Acts::Transform3::Identity();
    trafo.translate(Acts::Vector3(0., 0., positions[i] * cm));
    const auto rBounds = std::make_shared<const Acts::RadialBounds>(rmin, rmax);
    
    // creating custom surface array (surfaceArrayCreator.surfaceArrayOnDisc doesn't work)
    double tol = 1.;
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Closed> axisPhi(-M_PI, M_PI, nTubes);
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound> axisR(rmin, rmax, 1);
    using SGL = Acts::SurfaceArray::SurfaceGridLookup<decltype(axisR), decltype(axisPhi)>;
    std::vector<Acts::AxisDirection> axisDirections = {Acts::AxisDirection::AxisR, Acts::AxisDirection::AxisPhi};
    auto repsurface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rBounds);
    auto gridLookup = std::make_unique<SGL>(repsurface, tol, std::pair{axisR, axisPhi}, axisDirections);
    std::vector<const Acts::Surface*> surfacesRaw = unpackSmartPointers(vSurfaces);
    gridLookup->fill(gctx, surfacesRaw);
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(std::move(gridLookup),vSurfaces, trafo));

    auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), layerThickness, nullptr, Acts::active);

    for (auto& surface : vSurfaces) {
      auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
      mutableSurface->associateLayer(*layer);
    }

    layerOrderVec.emplace_back(layer, layer->referencePosition(gctx, Acts::AxisDirection::AxisZ));
  }
  auto binning = std::make_unique<const Acts::BinUtility>(layerBinEdges, Acts::open, Acts::AxisDirection::AxisZ, Acts::Transform3::Identity());
  auto layArr = std::make_unique<Acts::BinnedArrayXD<Acts::LayerPtr>>(layerOrderVec, std::move(binning));
 
  auto boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(2000._mm, 2000._mm, 3100._mm);
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), boundsVol, volumeMaterial, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, "Telescope");
  auto trackingGeometry = new Acts::TrackingGeometry(trackVolume);

  Acts::ObjVisualization3D objVis;
  const Acts::TrackingVolume& tgVolume = *(trackingGeometry->highestTrackingVolume());
  Acts::GeometryView3D::drawTrackingVolume(objVis, tgVolume, gctx);

  // Check geometry
  bool print_surface_info = 0;
  if (print_surface_info) {
    const Acts::TrackingVolume *highestTrackingVolume = trackingGeometry->highestTrackingVolume();
    printf("volumeId = %lu\n", highestTrackingVolume->geometryId().value());
    const Acts::LayerArray *confinedLayers = highestTrackingVolume->confinedLayers();
    for (const auto &layer : confinedLayers->arrayObjects())  {
      std::cout << "  layerId = " << layer->geometryId();
      printf("  thickness = %f type = %d\n", layer->thickness(), layer->layerType());
      if (layer->layerType()==-1) {
        const Acts::NavigationLayer* nlayer = dynamic_cast<const Acts::NavigationLayer*>(layer.get());
        //printf(" navigation %p\n", nlayer);
      } else if (layer->layerType()==1) {
        const Acts::DiscLayer* player = dynamic_cast<const Acts::DiscLayer*>(layer.get());
        printf(" disk %p\n", player);
        auto pos = Acts::Vector3(920,4,0);
        auto dir = Acts::Vector3(0,0,1);
        auto v = layer->surfaceArray()->neighbors(pos, dir);
        printf("    v.size()=%d\n",v.size());
        // for (auto& s : v) { std::cout << "    v.surfaceId = " << s->geometryId() << std::endl;  }
        // std::cout << std::endl;
      } else {
        const Acts::PlaneLayer* player = dynamic_cast<const Acts::PlaneLayer*>(layer.get());
        printf(" plane %p\n", player);
        auto pos = Acts::Vector3(920,4,0);
        auto dir = Acts::Vector3(0,0,1);
        auto v = layer->surfaceArray()->neighbors(pos, dir);
        printf("    v.size()=%d\n",v.size());
        // for (auto& s : v) { std::cout << "    v.surfaceId = " << s->geometryId() << std::endl;  }
        // std::cout << std::endl;
      }
      if (!layer->surfaceArray())
        continue;
      for (const auto &surface : layer->surfaceArray()->surfaces())
      {
        std::cout << "    surfaceId = " << surface->geometryId();
      }
    } //for layers
  } // print_surface_info


  return trackingGeometry;
}


Acts::TrackingGeometry* CreateTrackingGeometry2(bool addROC = 0, bool addFlange = 0){
  // Create materials
  // auto* nist = G4NistManager::Instance();
  // G4Material* siMat = nist->FindOrBuildMaterial("G4_Si");
  // G4Material* alMat = nist->FindOrBuildMaterial("G4_Al");
  // G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");

  // Acts::Geant4MaterialConverter converter;
  // Acts::Material silicon = converter.material(*siMat);
  // Acts::Material aluminium = converter.material(*alMat);
  // Acts::Material vacuum = converter.material(*worldMat);
  auto aluminium = Acts::Material::fromMassDensity(8.897_cm, 39.70_cm, 26.9815, 13, 2.699_g / 1_cm3);
  auto silicon = Acts::Material::fromMassDensity(9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  auto vacuum = Acts::Material::fromMassDensity(1.176e+4_cm, 7.204e+4_cm, 39.948, 18, 1.662e-10_g / 1_cm3);

  Acts::MaterialSlab matProp(silicon, thickness*cm);
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);
  const auto volumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(vacuum);

  // Construct the surfaces and layers
  Acts::LayerVector layVec;
 

  if (addROC) {
    double phi[nSectors];
    for (int i=0;i<nSectors;i++){
      phi[i] = (15.+30.*i)/180.*M_PI;
    }
    
    Acts::MaterialSlab matPropROC(aluminium, thicknessROC*cm);
    const auto surfaceMatROC = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropROC);
    Acts::Translation3 trans(0, 0, zROC*cm);
    Acts::Transform3 trafo(trans);
    const auto rBoundsROC = std::make_shared<const Acts::RadialBounds>(rMinROC*cm, rMaxROC*cm);
    auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMinROC*cm, rMaxROC*cm);
    surface->assignSurfaceMaterial(std::move(surfaceMatROC));
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
    auto layer = Acts::DiscLayer::create(trafo, rBoundsROC, std::move(surArray), thicknessROC*cm, nullptr, Acts::active);
    surface->associateLayer(*layer.get());
    layVec.push_back(layer);

    if (addFlange) {
      Acts::MaterialSlab matPropFrame(aluminium, thicknessFrame*cm);
      const auto surfaceMatFrame = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matPropFrame);
      std::vector<std::shared_ptr<const Acts::Surface>> vSurfaceFrame;

      Acts::Translation3 transFrameCircum1(0, 0, zFrameCircum1*cm);
      Acts::Translation3 transFrameCircum2(0, 0, zFrameCircum2*cm);
      Acts::Transform3 trafoFrameCircum1(transFrameCircum1);
      Acts::Transform3 trafoFrameCircum2(transFrameCircum2);

      const auto pBoundsFrameRadial = std::make_shared<const Acts::RectangleBounds>(halfXFrameRadial*cm, halfYFrameRadial*cm);
      for (int iSector = 0; iSector < nSectors; iSector++) {
        Acts::Transform3 trafoFrameRadial;
        trafoFrameRadial.setIdentity();
        trafoFrameRadial.rotate(Eigen::AngleAxisd(phi[iSector], Acts::Vector3(0, 0, 1)));
        trafoFrameRadial.translate(Acts::Vector3(cXFrameRadial0*cm, 0, zFrameRadial*cm));
        auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafoFrameRadial, pBoundsFrameRadial);
        surface->assignSurfaceMaterial(std::move(surfaceMatFrame));
        vSurfaceFrame.push_back(surface);

        auto bounds1 = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum1*cm, rMaxFrameCircum1*cm, M_PI/12.-eps, Acts::detail::radian_sym(iSector*M_PI/6));
        auto bounds2 = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum2*cm, rMaxFrameCircum2*cm, M_PI/12.-eps, Acts::detail::radian_sym(iSector*M_PI/6));
        auto surfaceFrameCircum1 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum1, bounds1);
        auto surfaceFrameCircum2 = Acts::Surface::makeShared<Acts::DiscSurface>(trafoFrameCircum2, bounds2);
        surfaceFrameCircum1->assignSurfaceMaterial(std::move(surfaceMatFrame));
        surfaceFrameCircum2->assignSurfaceMaterial(std::move(surfaceMatFrame));
        vSurfaceFrame.push_back(surfaceFrameCircum1);
        vSurfaceFrame.push_back(surfaceFrameCircum2);
      }

      Acts::SurfaceArrayCreator::Config sacConfig;
      Acts::SurfaceArrayCreator surfaceArrayCreator(sacConfig, Acts::getDefaultLogger("SurfaceArrayCreator", Acts::Logging::INFO));
      const auto rBoundsFrame = std::make_shared<const Acts::RadialBounds>(rMinFrameCircum1*cm, rMaxFrameCircum2*cm);
      Acts::Translation3 transFrame(0., 0., zFrameRadial*cm);
      Acts::Transform3 trafoFrame(transFrame);
      auto layerFrame = Acts::DiscLayer::create(trafoFrame, rBoundsFrame, std::move(surfaceArrayCreator.surfaceArrayOnDisc(gctx, vSurfaceFrame, 3, 12)), thicknessFrame*cm, nullptr, Acts::active);

      for (auto& surface : vSurfaceFrame) {
        auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
        mutableSurface->associateLayer(*layerFrame.get());
      }
      layVec.push_back(layerFrame);
    }
  }


  // const auto rBounds = std::make_shared<const Acts::RadialBounds>(rMinStation, rMaxStation);  // <- for disk-like layers
  const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rMaxStation*cm, rMaxStation*cm); // <- for square-like layers
  for (unsigned int i = 0; i < positions.size(); i++) {
    Acts::Transform3 trafo = Acts::Transform3::Identity();
    trafo.translate(Acts::Vector3(0, 0, positions[i]*cm));
    if (i%2==1) trafo.rotate(Eigen::AngleAxisd(std::numbers::pi/2, Acts::Vector3(0, 0, 1))); // for strip-like simulations
    // create surface
    // auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rMinStation, rMaxStation); // <- for disk-like layers
    auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds); // <- for square-like layers
    surface->assignSurfaceMaterial(std::move(surfaceMaterial));
    auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
    // create layer
    // auto layer = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), 1._mm); // <- for disk-like layers
    auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), 1._mm); // <- for square-like layers
    surface->associateLayer(*layer.get());
    layVec.push_back(layer);
    // create detector element
    auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thickness*cm);
    surface->assignDetectorElement(*detElement.get());
    detectorStore.push_back(std::move(detElement));
  }

  // Create layer array
  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(gctx, layVec, 1500_mm, positions.back()*cm + 2._mm, Acts::BinningType::arbitrary, Acts::AxisDirection::AxisZ));

  // Build mother tracking volume
  Acts::Translation3 transVol(0, 0, 0);
  Acts::Transform3 trafoVol(transVol);
  auto boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(2000._mm, 2000._mm, 3100._mm); // <- for square-like layers
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(trafoVol, boundsVol, volumeMaterial, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, "Telescope");
  // Build tracking geometry
  auto trackingGeometry = new Acts::TrackingGeometry(trackVolume);
  auto gctx = Acts::GeometryContext();
  Acts::ObjVisualization3D objVis;
  const Acts::TrackingVolume& tgVolume = *(trackingGeometry->highestTrackingVolume());
  Acts::GeometryView3D::drawTrackingVolume(objVis, tgVolume, gctx);


  // Check geometry
  bool print_surface_info = 0;
  if (print_surface_info) {
    const Acts::TrackingVolume *highestTrackingVolume = trackingGeometry->highestTrackingVolume();
    printf("volumeId = %lu\n", highestTrackingVolume->geometryId().value());
    const Acts::LayerArray *confinedLayers = highestTrackingVolume->confinedLayers();
    for (const auto &layer : confinedLayers->arrayObjects())  {
      printf("  layerId = %lu, thickness = %f type = %d\n", layer->geometryId().value(), layer->thickness(), layer->layerType());
      if (layer->layerType()==-1) {
        const Acts::NavigationLayer* nlayer = dynamic_cast<const Acts::NavigationLayer*>(layer.get());
        printf(" navigation %p\n", nlayer);
      } else if (layer->layerType()==0) {
        const Acts::DiscLayer* player = dynamic_cast<const Acts::DiscLayer*>(layer.get());
        printf(" disk %p\n", player);
      } else {
        const Acts::PlaneLayer* player = dynamic_cast<const Acts::PlaneLayer*>(layer.get());
        printf(" plane %p\n", player);
      }
      if (!layer->surfaceArray())
        continue;
      for (const auto &surface : layer->surfaceArray()->surfaces())
      {
        printf("    surfaceId = %lu\n", surface->geometryId().value());
      }
    } //for layers
  } // print_surface_info
  return trackingGeometry;
}

#endif