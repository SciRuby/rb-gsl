//TODO: More efficient sorting algorithm than selection sort would be nice, probably.
// Remember, we're dealing with unique keys, which simplifies things.
// Doesn't have to be in-place, since we probably just multiplied and that wasn't in-place.

void %%INT_ABBREV%%_%%REAL_ABBREV%%_smmp_sort_columns_(y_size_t n, YALE_PARAM A)
{
  u_%%INT%% i, jj, jj_start, jj_end, min, min_jj;
  %%REAL%% temp_val;

  u_%%INT%% *ia = (u_%%INT%%*)(A.ia),
            *ja = (u_%%INT%%*)(A.ja);
  %%REAL%% *a = (%%REAL%%*)(A.a);

  for (i = 0; i < n; ++i) {
    // No need to sort if there are 0 or 1 entries in the row
    if (ia[i+1] - ia[i] < 2) continue;

    jj_end = ia[i+1];
    for (jj_start = ia[i]; jj_start < jj_end; ++jj_start) {

      // If the previous min is just current-1, this key/value pair is already in sorted order.
      // This follows from the unique condition on our column keys.
      if (jj_start > ia[i] && min+1 == ja[jj_start]) {
        min    = ja[jj_start];
        continue;
      }

      // find the minimum key (column index) between jj_start and jj_end
      min    = ja[jj_start];
      min_jj = jj_start;
      for (jj = jj_start+1; jj < jj_end; ++jj) {
        if (ja[jj] < min) {
          min_jj = jj;
          min    = ja[jj];
        }
      }

      // if min is already first, skip this iteration
      if (min_jj == jj_start) continue;

      for (jj = jj_start; jj < jj_end; ++jj) {
        // swap minimum key/value pair with key/value pair in the first position.
        if (min_jj != jj) {
          // min already = ja[min_jj], so use this as temp_key
          temp_val = a[min_jj];

          ja[min_jj] = ja[jj];
           a[min_jj] =  a[jj];

          ja[jj] = min;
           a[jj] = temp_val;
        }
      }
    }
  }
}