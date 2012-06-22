inline TYPE div4(FLOAT xr, FLOAT xi, FLOAT yr, FLOAT yi) {
  double denom = yr * yr + xi * xi;
  return (struct TYPE) {
    (xr * yr + xi * yi) / denom,
    (xr * yi - xi * yr) / denom
  };
}
