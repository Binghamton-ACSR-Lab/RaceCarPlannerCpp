//
// Created by acsr on 11/22/22.
//

#ifndef RACECARPLANNER_DB_MANAGER_H
#define RACECARPLANNER_DB_MANAGER_H

#include "SQLiteCpp/SQLiteCpp.h"
#include <queue>
#include <thread>
#include <atomic>
#include <mutex>

using namespace std::chrono_literals;
namespace acsr {
    class DbManager {
    public:
        DbManager(const std::string &db_file) : db_(db_file, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE) {

            t = std::thread([this]() {
                while (run_flag) {
                    if (queue_.empty()) {
                        std::this_thread::sleep_for(200ms);
                        continue;
                    }
                    std::scoped_lock<std::mutex> lock(mutex);
                    db_.exec(queue_.front());
                    queue_.pop();
                }
            });

        }

        virtual ~DbManager() {
            run_flag = false;
        }

        void add_sql(const std::string &state) {
            std::scoped_lock<std::mutex> lock(mutex);
            queue_.push(state);
        }

        template<class T>
        void add_sql(T begin, T end) {
            std::scoped_lock lock(mutex);
            auto it = begin;
            while (it != end)
                queue_.push(*it++);
        }


    private:
        std::thread t;
        SQLite::Database db_;
        std::queue<std::string> queue_;
        std::atomic_bool run_flag = true;
        std::mutex mutex;

    };
}

#endif //RACECARPLANNER_DB_MANAGER_H
