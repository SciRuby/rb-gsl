inline TYPE mul4(INT anum, INT aden, INT bnum, INT bden) {
  INT g1 = gcf(anum, bden);
  INT g2 = gcf(aden, bnum);
  return (struct TYPE) { (anum / g1) * (bnum / g2), (aden / g2) * (bden / g1) };
}
