#ifndef HOSC_DEBUG_LOGGER_H
#define HOSC_DEBUG_LOGGER_H

#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <direct.h>

inline const char* hosc_debug_log_dir()
{
    return "D:\\Development\\DLL\\WithSGeMS\\Template\\log";
}

inline const char* hosc_debug_log_file()
{
    return "D:\\Development\\DLL\\WithSGeMS\\Template\\log\\HOSC_AVL_debug.log";
}

inline void hosc_debug_log(const char* stage, const char* fmt, ...)
{
    _mkdir(hosc_debug_log_dir());

    FILE* file = fopen(hosc_debug_log_file(), "a");
    if(file == NULL)
    {
        return;
    }

    time_t now = time(NULL);
    struct tm* local = localtime(&now);
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    int second = 0;
    if(local != NULL)
    {
        year = local->tm_year + 1900;
        month = local->tm_mon + 1;
        day = local->tm_mday;
        hour = local->tm_hour;
        minute = local->tm_min;
        second = local->tm_sec;
    }

    fprintf(
        file,
        "%04d-%02d-%02d %02d:%02d:%02d %s",
        year,
        month,
        day,
        hour,
        minute,
        second,
        stage == NULL ? "(null)" : stage);

    if(fmt != NULL && fmt[0] != '\0')
    {
        fprintf(file, " | ");
        va_list args;
        va_start(args, fmt);
        vfprintf(file, fmt, args);
        va_end(args);
    }

    fprintf(file, "\n");
    fflush(file);
    fclose(file);
}

inline void hosc_install_exception_filter()
{
    static bool installed = false;
    if(!installed)
    {
        installed = true;
        hosc_debug_log("logger.install", "debug logger initialized");
    }
}

#endif
