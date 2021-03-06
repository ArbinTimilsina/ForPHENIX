//
// File generated by rootcint at Tue Dec 31 08:18:43 2013

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME EMCalMap_Dict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "EMCalMap_Dict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
      #if !(defined(R__ACCESS_IN_SYMBOL) || defined(R__USE_SHADOW_CLASS))
      typedef ::EMCalMap EMCalMap;
      #else
      class EMCalMap  :  public ::SubsysReco {
         public:
         //friend XX;
         // To force the creation of a virtual table, throw just in case.
         virtual ~EMCalMap() throw() {};
          int verbo; //
          string outfname; //
         unsigned int runNumber; //
         unsigned int nRunEvents; //
         ::TH3D* hEmc3D[8]; //
         ::TFile* outfile; //
      };
      #endif

   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void EMCalMap_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void EMCalMap_Dictionary();
   static void delete_EMCalMap(void *p);
   static void deleteArray_EMCalMap(void *p);
   static void destruct_EMCalMap(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EMCalMap*)
   {
      // Make sure the shadow class has the right sizeof
      R__ASSERT(sizeof(::EMCalMap) == sizeof(::ROOT::Shadow::EMCalMap));
      ::EMCalMap *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::EMCalMap),0);
      static ::ROOT::TGenericClassInfo 
         instance("EMCalMap", "/direct/phenix+u/arbint/WarnMaps/code/CuAuEMCal/EMCalMap.h", 24,
                  typeid(::EMCalMap), DefineBehavior(ptr, ptr),
                  &EMCalMap_ShowMembers, &EMCalMap_Dictionary, isa_proxy, 4,
                  sizeof(::EMCalMap) );
      instance.SetDelete(&delete_EMCalMap);
      instance.SetDeleteArray(&deleteArray_EMCalMap);
      instance.SetDestructor(&destruct_EMCalMap);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EMCalMap*)
   {
      return GenerateInitInstanceLocal((::EMCalMap*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::EMCalMap*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void EMCalMap_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::EMCalMap*)0x0)->GetClass();
   }

} // end of namespace ROOT

//______________________________________________________________________________
namespace ROOT {
   void EMCalMap_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
      // Inspect the data members of an object of class EMCalMap.
      typedef ::ROOT::Shadow::EMCalMap ShadowClass;
      ShadowClass *sobj = (ShadowClass*)obj;
      if (sobj) { } // Dummy usage just in case there is no datamember.

      TClass *R__cl  = ::ROOT::GenerateInitInstanceLocal((const ::EMCalMap*)0x0)->GetClass();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "verbo", &sobj->verbo);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "outfname", (void*)&sobj->outfname);
      R__insp.InspectMember("const string", (void*)&sobj->outfname, "outfname.", false);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "runNumber", &sobj->runNumber);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "nRunEvents", &sobj->nRunEvents);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*hEmc3D[8]", &sobj->hEmc3D);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*outfile", &sobj->outfile);
      R__insp.GenericShowMembers("SubsysReco", ( ::SubsysReco * )( (::EMCalMap*) obj ), false);
   }

}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_EMCalMap(void *p) {
      delete ((::EMCalMap*)p);
   }
   static void deleteArray_EMCalMap(void *p) {
      delete [] ((::EMCalMap*)p);
   }
   static void destruct_EMCalMap(void *p) {
      typedef ::EMCalMap current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EMCalMap

/********************************************************
* EMCalMap_Dict.C
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableEMCalMap_Dict();

extern "C" void G__set_cpp_environmentEMCalMap_Dict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("/direct/phenix+u/arbint/WarnMaps/code/CuAuEMCal/EMCalMap.h");
  G__cpp_reset_tagtableEMCalMap_Dict();
}
#include <new>
extern "C" int G__cpp_dllrevEMCalMap_Dict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* EMCalMap */
static int G__EMCalMap_Dict_411_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   EMCalMap* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new EMCalMap(*((string*) G__int(libp->para[0])));
   } else {
     p = new((void*) gvp) EMCalMap(*((string*) G__int(libp->para[0])));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__EMCalMap_Dict_411_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((EMCalMap*) G__getstructoffset())->IsPbGl((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__EMCalMap_Dict_411_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   EMCalMap* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new EMCalMap(*(EMCalMap*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef EMCalMap G__TEMCalMap;
static int G__EMCalMap_Dict_411_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (EMCalMap*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((EMCalMap*) (soff+(sizeof(EMCalMap)*i)))->~G__TEMCalMap();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (EMCalMap*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((EMCalMap*) (soff))->~G__TEMCalMap();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* EMCalMap */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncEMCalMap_Dict {
 public:
  G__Sizep2memfuncEMCalMap_Dict(): p(&G__Sizep2memfuncEMCalMap_Dict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncEMCalMap_Dict::*p)();
};

size_t G__get_sizep2memfuncEMCalMap_Dict()
{
  G__Sizep2memfuncEMCalMap_Dict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceEMCalMap_Dict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap))) {
     EMCalMap *G__Lderived;
     G__Lderived=(EMCalMap*)0x1000;
     {
       SubsysReco *G__Lpbase=(SubsysReco*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap),G__get_linked_tagnum(&G__EMCalMap_DictLN_SubsysReco),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       Fun4AllBase *G__Lpbase=(Fun4AllBase*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap),G__get_linked_tagnum(&G__EMCalMap_DictLN_Fun4AllBase),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableEMCalMap_Dict() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__EMCalMap_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__EMCalMap_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__EMCalMap_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__EMCalMap_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<const std::string,TNamed*>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_maplEconstsPstringcOTNamedmUcOlesslEconstsPstringgRcOallocatorlEpairlEconstsPstringcOTNamedmUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<const string,TNamed*>",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_maplEconstsPstringcOTNamedmUcOlesslEconstsPstringgRcOallocatorlEpairlEconstsPstringcOTNamedmUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<const string,TNamed*,less<const string> >",117,G__get_linked_tagnum(&G__EMCalMap_DictLN_maplEconstsPstringcOTNamedmUcOlesslEconstsPstringgRcOallocatorlEpairlEconstsPstringcOTNamedmUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* EMCalMap */
static void G__setup_memvarEMCalMap(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap));
   { EMCalMap *p; p=(EMCalMap*)0x1000; if (p) { }
   G__memvar_setup((void*)0,105,0,1,-1,-1,-1,2,"verbo=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,1,G__get_linked_tagnum(&G__EMCalMap_DictLN_string),-1,-1,2,"outfname=",0,(char*)NULL);
   G__memvar_setup((void*)0,104,0,0,-1,-1,-1,2,"runNumber=",0,(char*)NULL);
   G__memvar_setup((void*)0,104,0,0,-1,-1,-1,2,"nRunEvents=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__EMCalMap_DictLN_TH3D),-1,-1,2,"hEmc3D[8]=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__EMCalMap_DictLN_TFile),-1,-1,2,"outfile=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarEMCalMap_Dict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncEMCalMap(void) {
   /* EMCalMap */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap));
   G__memfunc_setup("EMCalMap",704,G__EMCalMap_Dict_411_0_1, 105, G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap), -1, 0, 1, 1, 1, 0, "u 'string' - 0 - _outfilename", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Init",404,(G__InterfaceMethod) NULL,105, -1, -1, 0, 1, 1, 1, 0, "U 'PHCompositeNode' - 0 - topNode", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("InitRun",713,(G__InterfaceMethod) NULL,105, -1, -1, 0, 1, 1, 1, 0, "U 'PHCompositeNode' - 0 - topNode", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("process_event",1408,(G__InterfaceMethod) NULL,105, -1, -1, 0, 1, 1, 1, 0, "U 'PHCompositeNode' - 0 - topNode", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("End",279,(G__InterfaceMethod) NULL,105, -1, -1, 0, 1, 1, 1, 0, "U 'PHCompositeNode' - 0 - topNode", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("EndRun",588,(G__InterfaceMethod) NULL,105, -1, -1, 0, 1, 1, 1, 0, "i - - 10 - runNumber", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("IsPbGl",545,G__EMCalMap_Dict_411_0_7, 105, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - armsect", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("EMCalMap", 704, G__EMCalMap_Dict_411_0_8, (int) ('i'), G__get_linked_tagnum(&G__EMCalMap_DictLN_EMCalMap), -1, 0, 1, 1, 1, 0, "u 'EMCalMap' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~EMCalMap", 830, G__EMCalMap_Dict_411_0_9, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncEMCalMap_Dict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalEMCalMap_Dict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcEMCalMap_Dict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__EMCalMap_DictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_Fun4AllBase = { "Fun4AllBase" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_PHCompositeNode = { "PHCompositeNode" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_SubsysReco = { "SubsysReco" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_maplEconstsPstringcOTNamedmUcOlesslEconstsPstringgRcOallocatorlEpairlEconstsPstringcOTNamedmUgRsPgRsPgR = { "map<const string,TNamed*,less<const string>,allocator<pair<const string,TNamed*> > >" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_TH3D = { "TH3D" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_TFile = { "TFile" , 99 , -1 };
G__linked_taginfo G__EMCalMap_DictLN_EMCalMap = { "EMCalMap" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableEMCalMap_Dict() {
  G__EMCalMap_DictLN_string.tagnum = -1 ;
  G__EMCalMap_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__EMCalMap_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__EMCalMap_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__EMCalMap_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__EMCalMap_DictLN_Fun4AllBase.tagnum = -1 ;
  G__EMCalMap_DictLN_PHCompositeNode.tagnum = -1 ;
  G__EMCalMap_DictLN_SubsysReco.tagnum = -1 ;
  G__EMCalMap_DictLN_maplEconstsPstringcOTNamedmUcOlesslEconstsPstringgRcOallocatorlEpairlEconstsPstringcOTNamedmUgRsPgRsPgR.tagnum = -1 ;
  G__EMCalMap_DictLN_TH3D.tagnum = -1 ;
  G__EMCalMap_DictLN_TFile.tagnum = -1 ;
  G__EMCalMap_DictLN_EMCalMap.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableEMCalMap_Dict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_string);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_Fun4AllBase);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_PHCompositeNode);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_SubsysReco);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_maplEconstsPstringcOTNamedmUcOlesslEconstsPstringgRcOallocatorlEpairlEconstsPstringcOTNamedmUgRsPgRsPgR);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_TH3D);
   G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_TFile);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__EMCalMap_DictLN_EMCalMap),sizeof(EMCalMap),-1,295936,(char*)NULL,G__setup_memvarEMCalMap,G__setup_memfuncEMCalMap);
}
extern "C" void G__cpp_setupEMCalMap_Dict(void) {
  G__check_setup_version(30051515,"G__cpp_setupEMCalMap_Dict()");
  G__set_cpp_environmentEMCalMap_Dict();
  G__cpp_setup_tagtableEMCalMap_Dict();

  G__cpp_setup_inheritanceEMCalMap_Dict();

  G__cpp_setup_typetableEMCalMap_Dict();

  G__cpp_setup_memvarEMCalMap_Dict();

  G__cpp_setup_memfuncEMCalMap_Dict();
  G__cpp_setup_globalEMCalMap_Dict();
  G__cpp_setup_funcEMCalMap_Dict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncEMCalMap_Dict();
  return;
}
class G__cpp_setup_initEMCalMap_Dict {
  public:
    G__cpp_setup_initEMCalMap_Dict() { G__add_setup_func("EMCalMap_Dict",(G__incsetup)(&G__cpp_setupEMCalMap_Dict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initEMCalMap_Dict() { G__remove_setup_func("EMCalMap_Dict"); }
};
G__cpp_setup_initEMCalMap_Dict G__cpp_setup_initializerEMCalMap_Dict;

