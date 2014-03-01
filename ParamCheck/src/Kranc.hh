
class KrancData
{
public:
  int dir;
  int face;
  int imin[3];
  int imax[3];
  int tile_imin[3];
  int tile_imax[3];
  CCTK_REAL normal[3];
  CCTK_REAL tangentA[3];
  CCTK_REAL tangentB[3];
};

void ParamCheck_TiledLoopOverInterior(cGH const * restrict const cctkGH,
                                        void (calc)(const cGH* restrict const cctkGH, const KrancData &kd));

#define GFOffset(u,di,dj,dk) KRANC_GFOFFSET3D(&u[index],di,dj,dk)