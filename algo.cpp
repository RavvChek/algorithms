// Задача яндекс контеста по предмету алгоритмы и структуры данных в ИТМО №1

#include <iostream>
#include <vector>
using namespace std;

int main()
{
    int count;
    cin >> count;
    vector<int> flowers(count);
    for (int i = 0; i < count; i++) {
        cin >> flowers[i];
    }

    int prevNum = -1;
    int prevPrevNum = -1;
    int maxLength = 0;
    int curLength = 0;
    int index = 0;
    int startIndex = 1;
    int maxStartIndex = 1;
    for (int flower : flowers) {
        if (flower != prevNum || flower != prevPrevNum) {
            curLength++;
        }
        else {
            if (maxLength < curLength) {
                maxLength = curLength;
                maxStartIndex = startIndex;
            }
            curLength = 2;
            startIndex = index;
        }
        prevPrevNum = prevNum;
        prevNum = flower;
        index++;
    }
    if (maxLength < curLength) {
        maxLength = curLength;
        maxStartIndex = startIndex;
    }
    /*cout << maxLength << endl;*/
    cout << maxStartIndex << " ";
    cout << maxStartIndex + maxLength - 1;
}

// Задача яндекс контеста по предмету алгоритмы и структуры данных в ИТМО №2

