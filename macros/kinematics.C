void kinematics(){
  double r6 = 120;
  double eta6 = 1.14;
  double theta6 = 2 * TMath::ATan(TMath::Exp(-eta6));
  double tt6 = TMath::Tan(theta6);
  double z6 = 170;
//  double r = z6*tt6;
  double z = r6/tt6;
  printf("%f\n",z);


  return;

  double eta5 = 2.2;
  double theta5 = 2 * TMath::ATan(TMath::Exp(-eta5));
  double tt5 = TMath::Tan(theta5);
  double r5 = 2100*tt5;
  printf("%f\n", r5);
  

  double p4 = 2.35;
  double eta4 = 2.0;
  double theta4 = 2 * TMath::ATan(TMath::Exp(-eta4));
  double pt4 = p4*TMath::Sin(theta4);
  printf("%f\n", pt4);
  return;

  double eta = 1.6;
  TVector3 v;
  v.SetPtEtaPhi(0.35,eta,0);
  printf("%f\n", tan(v.Theta()));
  printf("%f\n", v.Mag());
  printf("%f\n", 1./v.Mag());
  return;


  double tanTheta = 1300/3000.;
  //double tanTheta = 340/1630.;
  TVector3 v2;
  v2.SetMagThetaPhi(0.2,atan(tanTheta),0);
  printf("%f\n", v2.Eta());

  double p3 = 0.4;
  double eta3 = 1.6;
  double theta3 = 2 * TMath::ATan(TMath::Exp(-eta3));
  double pt3 = p3*TMath::Sin(theta3);
  printf("%f\n", pt3);
  printf("%f\n",1./TMath::Tan(theta3));

}

