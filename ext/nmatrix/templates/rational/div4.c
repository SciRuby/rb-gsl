inline TYPE div(INT anum, INT aden, INT bnum, INT bden) {

  if (bnum < 0) {
    anum = -anum;
    bnum = -bnum;
  }

  return mul4(anum, aden, bden, bnum);
}
