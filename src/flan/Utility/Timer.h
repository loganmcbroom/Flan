#pragma once

#include <chrono>

namespace flan {

/** \cond */
class Timer
{
public:
    void start()
        {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
        }
    
    void stop()
        {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
        }
    
    float elapsed_ms()
        {
        std::chrono::time_point<std::chrono::system_clock> end_time;
        
        if(m_bRunning)
            end_time = std::chrono::system_clock::now();
        else
            end_time = m_EndTime;
        
        return std::chrono::duration_cast<std::chrono::milliseconds>( end_time - m_StartTime ).count();
     }
    
    float elapsedSeconds()
    {
        return elapsed_ms() / 1000.0f;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};
/** \endcond */

}