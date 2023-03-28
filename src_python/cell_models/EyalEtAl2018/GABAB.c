/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__GABAB
#define _nrn_initial _nrn_initial__GABAB
#define nrn_cur _nrn_cur__GABAB
#define _nrn_current _nrn_current__GABAB
#define nrn_jacob _nrn_jacob__GABAB
#define nrn_state _nrn_state__GABAB
#define _net_receive _net_receive__GABAB 
#define state state__GABAB 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define tau_r_GABAB _p[0]
#define tau_d_GABAB _p[1]
#define Use _p[2]
#define e _p[3]
#define i _p[4]
#define i_GABAB _p[5]
#define g_GABAB _p[6]
#define A_GABAB _p[7]
#define B_GABAB _p[8]
#define factor_GABAB _p[9]
#define DA_GABAB _p[10]
#define DB_GABAB _p[11]
#define v _p[12]
#define _g _p[13]
#define _tsav _p[14]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[1];
#define _gth 0
#define fGABAB_GABAB _thread1data[0]
#define fGABAB _thread[_gth]._pval[0]
#define u0 u0_GABAB
 double u0 = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau_r_GABAB", "ms",
 "tau_d_GABAB", "ms",
 "Use", "1",
 "e", "mV",
 "i", "nA",
 "i_GABAB", "nA",
 "g_GABAB", "uS",
 0,0
};
 static double A_GABAB0 = 0;
 static double B_GABAB0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "u0_GABAB", &u0_GABAB,
 "fGABAB_GABAB", &fGABAB_GABAB,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABAB",
 "tau_r_GABAB",
 "tau_d_GABAB",
 "Use",
 "e",
 0,
 "i",
 "i_GABAB",
 "g_GABAB",
 0,
 "A_GABAB",
 "B_GABAB",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 15, _prop);
 	/*initialize range parameters*/
 	tau_r_GABAB = 30;
 	tau_d_GABAB = 80;
 	Use = 1;
 	e = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 15;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GABAB_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 2,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 15, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 2;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABAB C:/Users/bherr035/OneDrive - Florida International University/!!PhD_Research/Python-files/Theta_paper_sim_mpi_implementation/cell_models/EyalEtAl2018/GABAB.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "GABAB receptor";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   DA_GABAB = - A_GABAB / tau_r_GABAB ;
   DB_GABAB = - B_GABAB / tau_d_GABAB ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 DA_GABAB = DA_GABAB  / (1. - dt*( ( - 1.0 ) / tau_r_GABAB )) ;
 DB_GABAB = DB_GABAB  / (1. - dt*( ( - 1.0 ) / tau_d_GABAB )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    A_GABAB = A_GABAB + (1. - exp(dt*(( - 1.0 ) / tau_r_GABAB)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_r_GABAB ) - A_GABAB) ;
    B_GABAB = B_GABAB + (1. - exp(dt*(( - 1.0 ) / tau_d_GABAB)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_d_GABAB ) - B_GABAB) ;
   }
  return 0;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[1] = _args[0] ;
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A_GABAB;
    double __primary = (A_GABAB + _args[1] * factor_GABAB) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau_r_GABAB ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau_r_GABAB ) - __primary );
    A_GABAB += __primary;
  } else {
 A_GABAB = A_GABAB + _args[1] * factor_GABAB ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B_GABAB;
    double __primary = (B_GABAB + _args[1] * factor_GABAB) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau_d_GABAB ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau_d_GABAB ) - __primary );
    B_GABAB += __primary;
  } else {
 B_GABAB = B_GABAB + _args[1] * factor_GABAB ;
     }
 } }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(1, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  A_GABAB = A_GABAB0;
  B_GABAB = B_GABAB0;
 {
   double _ltp_GABAB ;
 A_GABAB = 0.0 ;
   B_GABAB = 0.0 ;
   _ltp_GABAB = ( tau_r_GABAB * tau_d_GABAB ) / ( tau_d_GABAB - tau_r_GABAB ) * log ( tau_d_GABAB / tau_r_GABAB ) ;
   factor_GABAB = - exp ( - _ltp_GABAB / tau_r_GABAB ) + exp ( - _ltp_GABAB / tau_d_GABAB ) ;
   factor_GABAB = 1.0 / factor_GABAB ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   fGABAB = 33.33 * ( 0.5 - 2.0 / ( 1.0 + exp ( ( v + 98.73 ) / 12.5 ) ) ) ;
   g_GABAB = ( B_GABAB - A_GABAB ) * fGABAB ;
   i_GABAB = g_GABAB ;
   i = i_GABAB ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   state(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A_GABAB) - _p;  _dlist1[0] = &(DA_GABAB) - _p;
 _slist1[1] = &(B_GABAB) - _p;  _dlist1[1] = &(DB_GABAB) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "GABAB.mod";
static const char* nmodl_file_text = 
  "TITLE GABAB receptor\n"
  "\n"
  "\n"
  "COMMENT\n"
  "GABAB receptor conductance using a dual-exponential profile\n"
  "presynaptic short-term plasticity based on Fuhrmann et al, 2002\n"
  "Implemented by Beatriz Herrera April 8 based on the implementation by\n"
  "Srikanth Ramaswamy, Blue Brain Project, March 2009 and Guy 2018.\n"
  "Some modificaations of Ramaswamy implementation were made based on Yang et al.\n"
  "2016 Nat Commun. 7, 1\n"
  "\n"
  "\n"
  "14 (2016).\n"
  "GUY: Removed  plasticity and depression\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "NEURON {\n"
  "\n"
  "        POINT_PROCESS GABAB\n"
  "        RANGE  tau_r_GABAB, tau_d_GABAB\n"
  "        RANGE Use\n"
  "        RANGE i,  i_GABAB,  g_GABAB, e, gmax\n"
  "        NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "\n"
  "	   tau_r_GABAB = 30   (ms) : dual-exponential conductance profile\n"
  "        tau_d_GABAB = 80     (ms) : IMPORTANT: tau_r < tau_d\n"
  "        Use = 1.0   (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values)\n"
  "\n"
  "        e = 0     (mV)  : GABAB and NMDA reversal potential\n"
  "    	:gmax = .001 (uS) :1nS weight conversion factor (from nS to uS)\n"
  "    	u0 = 0 :initial value of u, which is the running value of Use\n"
  "        fGABAB\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1\n"
  "for comparison with Pr to decide whether to activate the synapse or not\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "\n"
  "        v (mV)\n"
  "        i (nA)\n"
  "	i_GABAB (nA)\n"
  "	g_GABAB (uS)\n"
  "	factor_GABAB\n"
  "\n"
  "}\n"
  "\n"
  "STATE {\n"
  "\n"
  "\n"
  "	A_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAB\n"
  "    B_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAB\n"
  "}\n"
  "\n"
  "INITIAL{\n"
  "\n"
  "    LOCAL  tp_GABAB\n"
  "\n"
  "	A_GABAB = 0\n"
  "	B_GABAB = 0\n"
  "\n"
  "	tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) :time to peak of the conductance\n"
  "\n"
  "\n"
  "\n"
  "	factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) :GABAB Normalization factor - so that when t = tp_GABAB, gsyn = gpeak\n"
  "    factor_GABAB = 1/factor_GABAB\n"
  "\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "\n"
  "    SOLVE state METHOD cnexp\n"
  "    fGABAB = 33.33*(0.5-2/(1+exp((v + 98.73)/12.5))) :fGABAB kinetics - Shoemaker, Neurocomputing 2011\n"
  "    g_GABAB = (B_GABAB-A_GABAB)*fGABAB :compute time varying conductance as the difference of state variables B_GABAB and A_GABAB and fGABAB kinetics\n"
  "\n"
  "	i_GABAB = g_GABAB :compute the GABAB driving force based on the time varying conductance, membrane potential, and GABAB reversal\n"
  "	i =  i_GABAB\n"
  "}\n"
  "\n"
  "DERIVATIVE state{\n"
  "\n"
  "\n"
  "	A_GABAB' = -A_GABAB/tau_r_GABAB\n"
  "    B_GABAB' = -B_GABAB/tau_d_GABAB\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "NET_RECEIVE (weight, weight_GABAB){\n"
  "\n"
  "	weight_GABAB = weight\n"
  "\n"
  "\n"
  "\n"
  "	A_GABAB = A_GABAB + weight_GABAB*factor_GABAB\n"
  "    B_GABAB = B_GABAB + weight_GABAB*factor_GABAB\n"
  "\n"
  "\n"
  "}\n"
  ;
#endif
