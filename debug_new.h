#ifndef _DEBUG_NEW_H_
#define _DEBUG_NEW_H_

    #ifdef _DEBUG

        #undef new
        extern void _RegDebugNew( void );
        extern void* __cdecl operator new( size_t, const char*, int );
        extern void __cdecl operator delete( void*, const char*, int);
        #define new new(__FILE__, __LINE__)
        
        #define REG_DEBUG_NEW _RegDebugNew();

    #else

        #define REG_DEBUG_NEW

    #endif // _DEBUG

#endif // _DEBUG_NEW_H_

