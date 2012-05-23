inline TYPE mul4(FLOAT xr, FLOAT xi, FLOAT yr, FLOAT yi) {
  return (struct TYPE) { xr * yr - xi * yi, xr * yi - xi * yr };
}