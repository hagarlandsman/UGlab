// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LinkDef_rdict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "LinkDef.cpp"
#include "UGburden.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_UGburden(void *p = nullptr);
   static void *newArray_UGburden(Long_t size, void *p);
   static void delete_UGburden(void *p);
   static void deleteArray_UGburden(void *p);
   static void destruct_UGburden(void *p);
   static void streamer_UGburden(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::UGburden*)
   {
      ::UGburden *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::UGburden >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("UGburden", ::UGburden::Class_Version(), "UGburden.h", 36,
                  typeid(::UGburden), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::UGburden::Dictionary, isa_proxy, 16,
                  sizeof(::UGburden) );
      instance.SetNew(&new_UGburden);
      instance.SetNewArray(&newArray_UGburden);
      instance.SetDelete(&delete_UGburden);
      instance.SetDeleteArray(&deleteArray_UGburden);
      instance.SetDestructor(&destruct_UGburden);
      instance.SetStreamerFunc(&streamer_UGburden);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::UGburden*)
   {
      return GenerateInitInstanceLocal(static_cast<::UGburden*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::UGburden*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr UGburden::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *UGburden::Class_Name()
{
   return "UGburden";
}

//______________________________________________________________________________
const char *UGburden::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UGburden*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int UGburden::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::UGburden*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *UGburden::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UGburden*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *UGburden::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::UGburden*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void UGburden::Streamer(TBuffer &R__b)
{
   // Stream an object of class UGburden.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> refEast;
      R__b >> refNorth;
      R__b.StreamObject(&(cloud),typeid(cloud));
      filename.Streamer(R__b);
      R__b >> kd_tree;
      R__b >> n_cols;
      R__b.CheckByteCount(R__s, R__c, UGburden::IsA());
   } else {
      R__c = R__b.WriteVersion(UGburden::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << refEast;
      R__b << refNorth;
      R__b.StreamObject(&(cloud),typeid(cloud));
      filename.Streamer(R__b);
      R__b << kd_tree;
      R__b << n_cols;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_UGburden(void *p) {
      return  p ? new(p) ::UGburden : new ::UGburden;
   }
   static void *newArray_UGburden(Long_t nElements, void *p) {
      return p ? new(p) ::UGburden[nElements] : new ::UGburden[nElements];
   }
   // Wrapper around operator delete
   static void delete_UGburden(void *p) {
      delete (static_cast<::UGburden*>(p));
   }
   static void deleteArray_UGburden(void *p) {
      delete [] (static_cast<::UGburden*>(p));
   }
   static void destruct_UGburden(void *p) {
      typedef ::UGburden current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_UGburden(TBuffer &buf, void *obj) {
      ((::UGburden*)obj)->::UGburden::Streamer(buf);
   }
} // end of namespace ROOT for class ::UGburden

namespace {
  void TriggerDictionaryInitialization_LinkDef_rdict_Impl() {
    static const char* headers[] = {
"LinkDef.cpp",
"UGburden.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/hagar/Tools/root/install_dir/include/",
"/home/hagar/Workspace/UGlab/kochav_Hayarden_work/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LinkDef_rdict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$UGburden.h")))  UGburden;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LinkDef_rdict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "LinkDef.cpp"
#include "UGburden.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"UGburden", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LinkDef_rdict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LinkDef_rdict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LinkDef_rdict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LinkDef_rdict() {
  TriggerDictionaryInitialization_LinkDef_rdict_Impl();
}
