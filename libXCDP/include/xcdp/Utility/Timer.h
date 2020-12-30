#pragma once

namespace xcdp {

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
    
    float elapsedMilliseconds()
        {
        std::chrono::time_point<std::chrono::system_clock> endTime;
        
        if(m_bRunning)
            endTime = std::chrono::system_clock::now();
        else
            endTime = m_EndTime;
        
        return std::chrono::duration_cast<std::chrono::milliseconds>( endTime - m_StartTime ).count();
     }
    
    float elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0f;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};
/** \endcond */

}