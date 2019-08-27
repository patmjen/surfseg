#include "util.h"

size_t std::hash<std::pair<int, int>>::operator()(const std::pair<int, int>& key) const noexcept
{
    static std::hash<long long> hasher;
    long long tohash = static_cast<long long>(key.first);
    tohash |= (static_cast<long long>(key.second) << 32);
    return hasher(tohash);
}

template <size_t N>
size_t std::hash<std::array<int, N>>::operator()(const std::array<int, N>& key) const noexcept
{
    // Code from https://codereview.stackexchange.com/q/171999
    static std::hash<int> hasher;
    size_t result = 0;
    for (int k : key) {
        result = result * 31 + hasher(k);
    }
    return result;
}

void graphErrFunc(const char *msg)
{
	// To avoid crashing the program we throw an exception here since Graph will otherwise call exit(1)
	throw std::runtime_error(msg);
}
