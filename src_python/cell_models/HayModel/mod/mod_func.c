#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _AMPA_reg();
extern void _CaDynamics_E2_reg();
extern void _Ca_HVA_reg();
extern void _Ca_LVAst_reg();
extern void _GABAA_reg();
extern void _GABABsyn_reg();
extern void _Ih_reg();
extern void _Im_reg();
extern void _K_Pst_reg();
extern void _K_Tst_reg();
extern void _NMDA_reg();
extern void _NaTa_t_reg();
extern void _NaTs2_t_reg();
extern void _Nap_Et2_reg();
extern void _SK_E2_reg();
extern void _SKv3_1_reg();
extern void _epsp_reg();
extern void _stim_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," AMPA.mod");
fprintf(stderr," CaDynamics_E2.mod");
fprintf(stderr," Ca_HVA.mod");
fprintf(stderr," Ca_LVAst.mod");
fprintf(stderr," GABAA.mod");
fprintf(stderr," GABABsyn.mod");
fprintf(stderr," Ih.mod");
fprintf(stderr," Im.mod");
fprintf(stderr," K_Pst.mod");
fprintf(stderr," K_Tst.mod");
fprintf(stderr," NMDA.mod");
fprintf(stderr," NaTa_t.mod");
fprintf(stderr," NaTs2_t.mod");
fprintf(stderr," Nap_Et2.mod");
fprintf(stderr," SK_E2.mod");
fprintf(stderr," SKv3_1.mod");
fprintf(stderr," epsp.mod");
fprintf(stderr," stim.mod");
fprintf(stderr, "\n");
    }
_AMPA_reg();
_CaDynamics_E2_reg();
_Ca_HVA_reg();
_Ca_LVAst_reg();
_GABAA_reg();
_GABABsyn_reg();
_Ih_reg();
_Im_reg();
_K_Pst_reg();
_K_Tst_reg();
_NMDA_reg();
_NaTa_t_reg();
_NaTs2_t_reg();
_Nap_Et2_reg();
_SK_E2_reg();
_SKv3_1_reg();
_epsp_reg();
_stim_reg();
}
