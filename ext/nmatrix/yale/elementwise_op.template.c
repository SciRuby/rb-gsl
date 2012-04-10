
static inline void %%TYPE_ABBREV%%_ew_op_unary(enum NMatrix_Ops op, %%TYPE%%* result, %%TYPE%% left) {
  switch(op) {
  case '!':
    %%TYPE *result = !left%%
    break;
  case NM_OP_NEG:
    %%TYPE *result = -left%%
    break;
  case '~':
    %%TYPE *result = ~left%%
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized element-wise unary operator");
  }
}



static inline void %%TYPE_ABBREV%%_ew_op_binary(enum NMatrix_Ops op, %%TYPE%%* result, %%TYPE%% left, %%TYPE%% right) {
  switch(op) {
  case '+':
    %%TYPE *result = left + right%%
    break;
  case '-':
    %%TYPE *result = left - right%%
    break;
  case '*':
    %%TYPE *result = left * right%%
    break;
  case '/':
    %%TYPE *result = left / right%%
    break;
  case '%':
    %%TYPE *result = left % right%%
    break;
  case NM_OP_EQEQ:
    %%TYPE *result = left == right%%
    break;
  case NM_OP_NEQ:
    %%TYPE *result = left != right%%
    break;
  case '>':
    %%TYPE *result = left > right%%
    break;
  case '<':
    %%TYPE *result = left < right%%
    break;
  case NM_OP_GTE:
    %%TYPE *result = left >= right%%
    break;
  case NM_OP_LTE:
    %%TYPE *result = left <= right%%
    break;
  case '&':
    %%TYPE *result = left & right%%
    break;
  case '|':
    %%TYPE *result = left | right%%
    break;
  case '^':
    %%TYPE *result = left ^ right%%
    break;
  case NM_OP_LSH:
    %%TYPE *result = left << right%%
    break;
  case NM_OP_RSH:
    %%TYPE *result = left >> right%%
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized element-wise binary operator");
  }
}
