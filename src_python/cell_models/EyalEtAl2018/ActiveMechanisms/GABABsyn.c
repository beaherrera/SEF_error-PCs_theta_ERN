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
 
#define nrn_init _nrn_init__GABAB3
#define _nrn_initial _nrn_initial__GABAB3
#define nrn_cur _nrn_cur__GABAB3
#define _nrn_current _nrn_current__GABAB3
#define nrn_jacob _nrn_jacob__GABAB3
#define nrn_state _nrn_state__GABAB3
#define _net_receive _net_receive__GABAB3 
#define bindkin bindkin__GABAB3 
 
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
#define gmax _p[0]
#define i _p[1]
#define g _p[2]
#define R _p[3]
#define Ron _p[4]
#define Roff _p[5]
#define G _p[6]
#define Gn _p[7]
#define edc _p[8]
#define synon _p[9]
#define Rinf _p[10]
#define Rtau _p[11]
#define Beta _p[12]
#define DRon _p[13]
#define DRoff _p[14]
#define DG _p[15]
#define v _p[16]
#define _g _p[17]
#define _tsav _p[18]
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
#define Cdur Cdur_GABAB3
 double Cdur = 5;
#define Cmax Cmax_GABAB3
 double Cmax = 0.5;
#define Erev Erev_GABAB3
 double Erev = -95;
#define KD KD_GABAB3
 double KD = 100;
#define K4 K4_GABAB3
 double K4 = 0.033;
#define K3 K3_GABAB3
 double K3 = 0.098;
#define K2 K2_GABAB3
 double K2 = 0.0013;
#define K1 K1_GABAB3
 double K1 = 0.52;
#define n n_GABAB3
 double n = 4;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Cmax_GABAB3", "mM",
 "Cdur_GABAB3", "ms",
 "K1_GABAB3", "/ms",
 "K2_GABAB3", "/ms",
 "K3_GABAB3", "/ms",
 "K4_GABAB3", "/ms",
 "Erev_GABAB3", "mV",
 "gmax", "uS",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double G0 = 0;
 static double Roff0 = 0;
 static double Ron0 = 0;
 static double delta_t = 1;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Cmax_GABAB3", &Cmax_GABAB3,
 "Cdur_GABAB3", &Cdur_GABAB3,
 "K1_GABAB3", &K1_GABAB3,
 "K2_GABAB3", &K2_GABAB3,
 "K3_GABAB3", &K3_GABAB3,
 "K4_GABAB3", &K4_GABAB3,
 "KD_GABAB3", &KD_GABAB3,
 "n_GABAB3", &n_GABAB3,
 "Erev_GABAB3", &Erev_GABAB3,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABAB3",
 "gmax",
 0,
 "i",
 "g",
 "R",
 0,
 "Ron",
 "Roff",
 "G",
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
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	gmax = 0.0001;
  }
 	_prop->param = _p;
 	_prop->param_size = 19;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
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
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GABABsyn_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 19, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 3;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABAB3 C:/Users/bherr035/OneDrive - Florida International University/!!PhD_Research/Python-files/Theta_paper_sim_mpi_implementation/cell_models/EyalEtAl2018/ActiveMechanisms/GABABsyn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int bindkin(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   DRon = synon * K1 * Cmax - ( K1 * Cmax + K2 ) * Ron ;
   DRoff = - K2 * Roff ;
   R = Ron + Roff ;
   DG = K3 * R - K4 * G ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 DRon = DRon  / (1. - dt*( ( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ) )) ;
 DRoff = DRoff  / (1. - dt*( ( - K2 )*( 1.0 ) )) ;
 R = Ron + Roff ;
 DG = DG  / (1. - dt*( ( - ( K4 )*( 1.0 ) ) )) ;
  return 0;
}
 /*END CVODE*/
 static int bindkin (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    Ron = Ron + (1. - exp(dt*(( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ))))*(- ( ( ( synon )*( K1 ) )*( Cmax ) ) / ( ( - ( ( ( K1 )*( Cmax ) + K2 ) )*( 1.0 ) ) ) - Ron) ;
    Roff = Roff + (1. - exp(dt*(( - K2 )*( 1.0 ))))*(- ( 0.0 ) / ( ( - K2 )*( 1.0 ) ) - Roff) ;
   R = Ron + Roff ;
    G = G + (1. - exp(dt*(( - ( K4 )*( 1.0 ) ))))*(- ( ( K3 )*( R ) ) / ( ( - ( K4 )*( 1.0 ) ) ) - G) ;
   }
  return 0;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 1.0 ) {
     _args[1] = _args[0] * ( Rinf + ( _args[1] - Rinf ) * exp ( - ( t - _args[2] ) / Rtau ) ) ;
     _args[2] = t ;
     synon = synon - _args[0] ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron - _args[1] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ) ) ) )*( - ( ( ( synon )*( K1 ) )*( Cmax ) ) / ( ( - ( ( ( K1 )*( Cmax ) + K2 ) )*( 1.0 ) ) ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron - _args[1]  ;
       }
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff + _args[1] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - K2 )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - K2 )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff + _args[1]  ;
       }
 }
   else {
     _args[1] = _args[0] * _args[1] * exp ( - Beta * ( t - _args[2] ) ) ;
     _args[2] = t ;
     synon = synon + _args[0] ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Ron;
    double __primary = (Ron + _args[1] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - ( ( K1 * Cmax + K2 ) )*( 1.0 ) ) ) ) )*( - ( ( ( synon )*( K1 ) )*( Cmax ) ) / ( ( - ( ( ( K1 )*( Cmax ) + K2 ) )*( 1.0 ) ) ) - __primary );
    Ron += __primary;
  } else {
 Ron = Ron + _args[1]  ;
       }
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = Roff;
    double __primary = (Roff - _args[1] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - K2 )*( 1.0 ) ) ) )*( - ( 0.0 ) / ( ( - K2 )*( 1.0 ) ) - __primary );
    Roff += __primary;
  } else {
 Roff = Roff - _args[1]  ;
       }
 net_send ( _tqitem, _args, _pnt, t +  Cdur , 1.0 ) ;
     }
   } }
 
static int _ode_count(int _type){ return 3;}
 
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
	for (_i=0; _i < 3; ++_i) {
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

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  G = G0;
  Roff = Roff0;
  Ron = Ron0;
 {
   R = 0.0 ;
   G = 0.0 ;
   synon = 0.0 ;
   Rinf = K1 * Cmax / ( K1 * Cmax + K2 ) ;
   Rtau = 1.0 / ( K1 * Cmax + K2 ) ;
   Beta = K2 ;
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
   Gn = G * G * G * G ;
   g = gmax * Gn / ( Gn + KD ) ;
   i = g * ( v - Erev ) ;
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
 {   bindkin(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(Ron) - _p;  _dlist1[0] = &(DRon) - _p;
 _slist1[1] = &(Roff) - _p;  _dlist1[1] = &(DRoff) - _p;
 _slist1[2] = &(G) - _p;  _dlist1[2] = &(DG) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "GABABsyn.mod";
static const char* nmodl_file_text = 
  ": $Id: gabab.mod,v 1.9 2004/06/17 16:04:05 billl Exp $\n"
  "\n"
  "COMMENT\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "	Kinetic model of GABA-B receptors\n"
  "	=================================\n"
  "\n"
  "  MODEL OF SECOND-ORDER G-PROTEIN TRANSDUCTION AND FAST K+ OPENING\n"
  "  WITH COOPERATIVITY OF G-PROTEIN BINDING TO K+ CHANNEL\n"
  "\n"
  "  PULSE OF TRANSMITTER\n"
  "\n"
  "  SIMPLE KINETICS WITH NO DESENSITIZATION\n"
  "\n"
  "	Features:\n"
  "\n"
  "  	  - peak at 100 ms; time course fit to Tom Otis' PSC\n"
  "	  - SUMMATION (psc is much stronger with bursts)\n"
  "\n"
  "\n"
  "	Approximations:\n"
  "\n"
  "	  - single binding site on receptor\n"
  "	  - model of alpha G-protein activation (direct) of K+ channel\n"
  "	  - G-protein dynamics is second-order; simplified as follows:\n"
  "		- saturating receptor\n"
  "		- no desensitization\n"
  "		- Michaelis-Menten of receptor for G-protein production\n"
  "		- \"resting\" G-protein is in excess\n"
  "		- Quasi-stat of intermediate enzymatic forms\n"
  "	  - binding on K+ channel is fast\n"
  "\n"
  "\n"
  "	Kinetic Equations:\n"
  "\n"
  "	  dR/dt = K1 * T * (1-R-D) - K2 * R\n"
  "\n"
  "	  dG/dt = K3 * R - K4 * G\n"
  "\n"
  "	  R : activated receptor\n"
  "	  T : transmitter\n"
  "	  G : activated G-protein\n"
  "	  K1,K2,K3,K4 = kinetic rate cst\n"
  "\n"
  "  n activated G-protein bind to a K+ channel:\n"
  "\n"
  "	n G + C <-> O		(Alpha,Beta)\n"
  "\n"
  "  If the binding is fast, the fraction of open channels is given by:\n"
  "\n"
  "	O = G^n / ( G^n + KD )\n"
  "\n"
  "  where KD = Beta / Alpha is the dissociation constant\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "  Parameters estimated from patch clamp recordings of GABAB PSP's in\n"
  "  rat hippocampal slices (Otis et al, J. Physiol. 463: 391-407, 1993).\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "  PULSE MECHANISM\n"
  "\n"
  "  Kinetic synapse with release mechanism as a pulse.\n"
  "\n"
  "  Warning: for this mechanism to be equivalent to the model with diffusion\n"
  "  of transmitter, small pulses must be used...\n"
  "\n"
  "  For a detailed model of GABAB:\n"
  "\n"
  "  Destexhe, A. and Sejnowski, T.J.  G-protein activation kinetics and\n"
  "  spill-over of GABA may account for differences between inhibitory responses\n"
  "  in the hippocampus and thalamus.  Proc. Natl. Acad. Sci. USA  92:\n"
  "  9515-9519, 1995.\n"
  "\n"
  "  For a review of models of synaptic currents:\n"
  "\n"
  "  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of\n"
  "  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition;\n"
  "  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1996.\n"
  "\n"
  "  This simplified model was introduced in:\n"
  "\n"
  "  Destexhe, A., Bal, T., McCormick, D.A. and Sejnowski, T.J.\n"
  "  Ionic mechanisms underlying synchronized oscillations and propagating\n"
  "  waves in a model of ferret thalamic slices. Journal of Neurophysiology\n"
  "  76: 2049-2070, 1996.\n"
  "\n"
  "  See also http://www.cnl.salk.edu/~alain\n"
  "\n"
  "\n"
  "\n"
  "  Alain Destexhe, Salk Institute and Laval University, 1995\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS GABAB3\n"
  "	RANGE R, G, g, gmax\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL Cmax, Cdur\n"
  "	GLOBAL K1, K2, K3, K4, KD, Erev\n"
  "}\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "     gmax = 0.0001 (uS)\n"
  "	Cmax	= 0.5	(mM)		: max transmitter concentration\n"
  "	Cdur	= 5	(ms)		: transmitter duration (rising phase)\n"
  ":\n"
  ":	From Kfit with long pulse (5ms 0.5mM)\n"
  ":\n"
  "	K1	= 0.52	(/ms mM)	: forward binding rate to receptor\n"
  "	K2	= 0.0013 (/ms)		: backward (unbinding) rate of receptor\n"
  "	K3	= 0.098 (/ms)		: rate of G-protein production\n"
  "	K4	= 0.033 (/ms)		: rate of G-protein decay\n"
  "	KD	= 100			: dissociation constant of K+ channel\n"
  "	n	= 4			: nb of binding sites of G-protein on K+\n"
  "	Erev	= -95	(mV)		: reversal potential (E_K)\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)		: postsynaptic voltage\n"
  "	i 		(nA)		: current = g*(v - Erev)\n"
  "	g 		(umho)		: conductance\n"
  "	Gn\n"
  "	R				: fraction of activated receptor\n"
  "	edc\n"
  "	synon\n"
  "	Rinf\n"
  "	Rtau (ms)\n"
  "	Beta (/ms)\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "	Ron Roff\n"
  "	G				: fraction of activated G-protein\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL {\n"
  "	R = 0\n"
  "	G = 0\n"
  "	synon = 0\n"
  "	Rinf = K1*Cmax/(K1*Cmax + K2)\n"
  "	Rtau = 1/(K1*Cmax + K2)\n"
  "	Beta = K2\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE bindkin METHOD cnexp\n"
  "	Gn = G*G*G*G : ^n = 4\n"
  "	g = gmax * Gn / (Gn+KD)\n"
  "	i = g*(v - Erev)\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE bindkin {\n"
  "	Ron' = synon*K1*Cmax - (K1*Cmax + K2)*Ron\n"
  "	Roff' = -K2*Roff\n"
  "	R = Ron + Roff\n"
  "	G' = K3 * R - K4 * G\n"
  "}\n"
  "\n"
  ": following supports both saturation from single input and\n"
  ": summation from multiple inputs\n"
  ": Note: automatic initialization of all reference args to 0\n"
  ": except first\n"
  "\n"
  "NET_RECEIVE(weight,  r0, t0 (ms)) {\n"
  "	if (flag == 1) { : at end of Cdur pulse so turn off\n"
  "		r0 = weight*(Rinf + (r0 - Rinf)*exp(-(t - t0)/Rtau))\n"
  "		t0 = t\n"
  "		synon = synon - weight\n"
  "		state_discontinuity(Ron, Ron - r0)\n"
  "		state_discontinuity(Roff, Roff + r0)\n"
  "        }else{ : at beginning of Cdur pulse so turn on\n"
  "		r0 = weight*r0*exp(-Beta*(t - t0))\n"
  "		t0 = t\n"
  "		synon = synon + weight\n"
  "		state_discontinuity(Ron, Ron + r0)\n"
  "		state_discontinuity(Roff, Roff - r0)\n"
  "		:come again in Cdur\n"
  "		net_send(Cdur, 1)\n"
  "        }\n"
  "}\n"
  ;
#endif
