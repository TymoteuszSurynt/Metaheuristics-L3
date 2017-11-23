#include <iostream>
#include <vector>
#include <cmath>
#include <random>
//#include <fstream>
#include <ctime>
#include <algorithm>

using namespace std;
typedef mt19937 rng_type;


rng_type rng;
//struct
struct city {
    double x;//coordinates
    double y;//coordinates
    double *distance; //table of distances to any other city
    bool visited;
};
struct individual{
    double path;
    int *permutation;
};
//Utils
void insertionSort1(individual *array,int start, int end){
    int j;
    for(int i=start+1;i<=end;i++){
        j=i;
        while(j>start){
            if(array[j].path<array[j-1].path) {
                swap(array[j], array[j-1]);
                j--;
            }else{
                break;
            }
        }
    }
}
void mergeMS(individual *a, int start, int mid, int end){
    individual *array2=new individual[end+1];
    int istart=start,jstart=mid+1,kstart=start;

    while (istart <= mid && jstart <= end){
        if (a[istart].path < a[jstart].path){
            array2[kstart] = a[istart];
            kstart++;
            istart++;
        }else{
            array2[kstart] = a[jstart];
            kstart++;
            jstart++;
        }
    }
    while (istart <= mid){
        array2[kstart] = a[istart];
        kstart++;
        istart++;
    }
    while (jstart <= end){
        array2[kstart] = a[jstart];
        kstart++;
        jstart++;
    }
    for (int i=start; i< kstart;i++){
        a[i] = array2[i];
    }
    delete [] array2;
}
void mergeSortInsert(individual* array, int start, int end){
    if(end-start<8){
        insertionSort1(array,start,end);
    }else {
        int mid = (start + end) / 2;
        mergeSortInsert(array, start, mid);
        mergeSortInsert(array, mid + 1, end);
        mergeMS(array, start, mid, end);
    }
}
void mergeSortInsert(individual* array,int size){
    mergeSortInsert(array,0,size-1);
}
int partition(individual * array, int start, int end)
{
    int i=start-1;
    double pivot=array[end].path;
    for(int j=start;j<=end-1;j++){
        if(array[j].path<=pivot){
            i++;

            swap(array[i],array[j]);
        }
    }

    swap(array[i+1],array[end]);
    return i+1;
}

void quickSortInsert(individual* array, int start, int end) {
    if(end-start<8){
        insertionSort1(array,start,end);
    }else {
        if (start < end) {
            int p = partition(array, start, end);
            quickSortInsert(array, start, p - 1);
            quickSortInsert(array, p + 1, end);
        }
    }
}
//print path
void printGraph(int *result, int size) {
    int check = -1;
    for (int i = 0; i < size; i++) {
        if (check != -1) {
            cerr << result[i] + 1 << endl;
        } else {
            if (result[i] == 0) {
                check = i;
                cerr << result[i] + 1 << endl;
            }
        }
    }
    for (int j = 0; j < check; ++j) {
        cerr << result[j] + 1 << endl;
    }
    cerr << 1 << endl;
}
//sum of distances in the result
double sumOfDistances(int *result, vector<struct city> &cities, int size) {
    double sum = 0;
    for (int i = 0; i < size - 1; i++) {
        sum += cities.at((unsigned long long int) result[i]).distance[result[i + 1]];
    }
    sum += cities.at((unsigned long long int) result[size - 1]).distance[result[0]];
    return sum;
}
//Counting distance
double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));
}
//function allowing to change 2 nodes in the graph
void opt2(int *result, int *result2, int size, int x, int y) {
    int i, j = 0;
    if (result[x] != result[y]) {
        for (i = 0; i < x; i++) {
            result2[i] = result[i];
        }
        for (i = x; i < y + 1; i++) {
            result2[i] = result[y - j];
            j++;
        }
        for (i = y + 1; i < size; i++) {
            result2[i] = result[i];
        }
    }
}
void shuffle(int* permutation,int a) {
    int j;
    for(int i=0;i<a-2;++i){
        j=rand()%(a-i)+i;
        swap(permutation[i],permutation[j]);
    }
}
//generate permutation
struct individual generateIndividual(int a,vector<struct city> &cities,int* result){
    struct individual individual;
    individual.permutation=new int[a];
    uniform_int_distribution<int> udist(0, a-1);
    int num1,num2;
    //int *result2 = new int[a];*/
    //for (int i = 0; i < a; ++i) {
    //   individual.permutation[i]=result[i];
    //}
    /*for (int l = 0; l < log(a)*20; ++l) {
        num1=udist(rng);
        num2=udist(rng);
        while(num1==num2){
            num2=udist(rng);
        }
        opt2(individual.permutation,result2,a,num1,num2);
        for (int k = 0; k < a; ++k) {
            individual.permutation[k]=result2[k];
        }
    }*/
    //shuffle(individual.permutation,a);
    num1=udist(rng);
    num2=udist(rng);
    while(num1==num2){
        num2=udist(rng);
    }
    opt2(result,individual.permutation,a,num1,num2);
    individual.path=sumOfDistances(individual.permutation,cities,a);
    return individual;
}



//function for reproduction
individual getChild(individual p1, individual p2, int size, vector<struct city> &cities){
    int ** tab=new int * [size];
    for (int i = 0; i < size; ++i) {
        tab[i]=new int[4];
        tab[i][0]=-1;
        tab[i][1]=-1;
        tab[i][2]=-1;
        tab[i][3]=-1;
    }
    individual child;
    vector <int> list;
    child.permutation=new int[size];
    int num;
    tab[0][0]=p1.permutation[size-1];
    tab[0][1]=p1.permutation[1];
    tab[0][2]=p2.permutation[size-1];
    tab[0][3]=p2.permutation[1];

    tab[size-1][0]=p1.permutation[size-2];
    tab[size-1][1]=p1.permutation[0];
    tab[size-1][2]=p2.permutation[size-2];
    tab[size-1][3]=p2.permutation[0];
    list.push_back(0);
    list.push_back(size-1);
    for (int i = 1; i < size-1; ++i) {
        num=p1.permutation[i];
        if(tab[num][0]==-1){
            tab[num][0]=p1.permutation[i-1];
            tab[num][1]=p1.permutation[i+1];
        }else if(tab[num][2]==-1){
            tab[num][2]=p1.permutation[i-1];
            tab[num][3]=p1.permutation[i+1];
        }
        num=p2.permutation[i];
        if(tab[num][0]==-1){
            tab[num][0]=p1.permutation[i-1];
            tab[num][1]=p1.permutation[i+1];
        }else if(tab[num][2]==-1){
            tab[num][2]=p1.permutation[i-1];
            tab[num][3]=p1.permutation[i+1];
        }
        list.push_back(i);
    }
    bool check;
    int city,newCity=0,h=1;
    int r= rand() % size;
    city=list.at((unsigned long long int) r);
    child.permutation[0]=city;
    list.erase(list.begin()+r);
    while(h<size){
        for (int i = 0; i < size; ++i) {
            if(tab[i][0]==city){
                tab[i][0]=-1;
            }
            if(tab[i][1]==city){
                tab[i][1]=-1;
            }
            if(tab[i][2]==city){
                tab[i][2]=-1;
            }
            if(tab[i][3]==city){
                tab[i][3]=-1;
            }
        }
        check=true;
        for (int j = 0; j < 4; ++j) {
            if(tab[city][j] != -1) {
                for (int i = j; i < 4; ++i) {
                    if (tab[city][i] == tab[city][j]) {
                        newCity = tab[city][i];
                        check = false;
                        break;
                    }
                }
            }
        }
        if(check){
            int dist=0,candidate=-1;
            for (int i = 0; i < 4; ++i) {
                if(tab[city][i]!=-1){
                    if(candidate!=-1){
                        if(dist>cities.at((unsigned long long int) tab[city][i]).distance[child.permutation[h - 1]]){
                            dist= (int) cities.at(
                                    (unsigned long long int) tab[city][i]).distance[child.permutation[h - 1]];
                            candidate=tab[city][i];
                        }
                    }else{
                        dist= (int) cities.at((unsigned long long int) tab[city][i]).distance[child.permutation[h - 1]];
                        candidate=tab[city][i];
                    }
                }
            }
            if(candidate!=-1){
                city=candidate;
                child.permutation[h]=city;
                list.erase(remove(list.begin(), list.end(), city), list.end());
                h++;
            }else {
                if(h==size){
                    break;
                }else {
                    r= rand() % (size-h);
                    city=list.at((unsigned long long int) r);
                    list.erase(list.begin()+r);
                    child.permutation[h]=city;
                    h++;
                }
            }
        }else{
            city=newCity;
            child.permutation[h]=city;
            h++;
            list.erase(remove(list.begin(), list.end(), city), list.end());
        }

    }
    uniform_int_distribution<int> udist(0, size-1);
    uniform_int_distribution<int> udist1(0, (int) log(size));
    int num1, num2;//,num3=udist1(rng);
//    for (int l = 0; l < num3; ++l) {
//        num1=udist(rng);
//        num2=udist(rng);
//        while(num1==num2){
//            num2=udist(rng);
//        }
//        swap(child.permutation[num1],child.permutation[num2]);
//    }
    num1=udist(rng);
    num2=udist(rng);
    int *result=new int[size];
    while(num1==num2){
        num2=udist(rng);
    }
    opt2(child.permutation,result,size,num1,num2);
    for (int k = 0; k < size; ++k) {
        child.permutation[k]=result[k];
    }
    child.path=sumOfDistances(child.permutation,cities,size);
    return child;

}
//function for first path
int greedyDistance(struct city city, int size, vector<struct city> &cities) {
    int min = -1;
    unsigned i = 0;
    double minValue = -1;
    for (; i < size; i++) {
        if (city.distance[i] != 0) {
            if (!cities.at(i).visited) {
                minValue = city.distance[i];
                min = i;
                break;
            }
        }
    }
    if (minValue == -1) {
        return -1;
    }
    for (; i < size; ++i) {
        if (city.distance[i] < minValue && city.distance[i] != 0) {
            if (!cities.at(i).visited) {
                minValue = city.distance[i];
                min = i;
            }
        }
    }
    return min;
}
int main() {
    clock_t begin = clock();
//    std::ifstream in("input.txt");
//    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
//    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!
    srand((unsigned int) time(NULL));

    int a, cityName, timeSec, iterations, population, change;
    double xCity, yCity,prevSum;
    vector<struct city> cities;
    cin >> a;
    population= (int) (log(a));
    population *= 20;
    for (int i = 0; i < a; i++) {
        struct city city;
        cin >> cityName;
        cin >> xCity >> yCity;
        city.x = xCity;
        city.y = yCity;
        city.visited=false;
        cities.push_back(city);
    }
    cin >> timeSec;
    if(timeSec>40){
        timeSec-=30;
    }
    for (unsigned i = 0; i < a; i++) {
        cities.at(i).distance = new double[a];
        for (int j = 0; j < a; j++) {
            cities.at(i).distance[j] = distance(cities[i].x, cities[i].y, cities[j].x, cities[j].y);
        }
    }
    double pathSum=0;
    int b=0,w=1;
    int *result = new int[a];
    result[0]=0;
    while (true) {
        cities.at((unsigned long long int) b).visited = true;
        b= greedyDistance(cities.at((unsigned long long int) b), a, cities);
        if (b == -1) {
            pathSum += cities.at((unsigned long long int) result[w - 1]).distance[0];
            break;
        } else {
            result[w] = b;
            pathSum += cities.at((unsigned long long int) result[w - 1]).distance[result[w]];
            w++;
        }

    }
    individual* individuals=new individual[population];
    individuals[0].path=pathSum;
    individuals[0].permutation=result;
    iterations = (int) log(a) * 8;
    for (int k = 1; k < population; ++k) {
        individuals[k]=generateIndividual(a,cities,result);
    }
    int r,r1=0;
    uniform_int_distribution<int> udist(0, population/2);
    prevSum=0;
    change= (int) log(a);
    while (iterations > 0) {
        quickSortInsert(individuals,0,population-1);

        for (int i = population/2; i < population; ++i) {
            r=udist(rng);
            while(r==r1){
                r=udist(rng);
            }
            individuals[i]=getChild(individuals[r1],individuals[r],a,cities);
            r1=udist(rng);
        }
        if(prevSum==individuals[0].path){
            change--;
            if(a==0){
                break;
            }
        }else{
            change=(int) log(a)*40;
        }
        if ((double(clock() - begin) / CLOCKS_PER_SEC) > timeSec) {
            break;
        }
        --iterations;
//cout<<(double(clock() - begin) / CLOCKS_PER_SEC)<<endl;
    }
    quickSortInsert(individuals,0,population-1);
    cout << individuals[0].path <<endl;
    printGraph(individuals[0].permutation, a);
    for (unsigned i = 0; i < a; i++) {
        delete[] cities.at(i).distance;
    }
    delete[] individuals;
    return 0;
}