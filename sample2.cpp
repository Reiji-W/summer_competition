#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <cfloat>

#define CHROM_LENGTH        500     // 遺伝子の長さ、商品数と同じ
#define GENERATION_MAX		350	    // 世代交代数
#define POP_SIZE			170    // 個体集団内の個体の数
#define ELITE               1
#define PROBABILITY         0.001


//本番の評価関数でも下記と同じ値を使います
#define CAL_LIMIT_PER_DAY   2400      //一日分の食事の摂取カロリーの上限、一日の区切りを決めるのに使用、勘で決めているので根拠は無い
#define IDEAL_CAL_PER_DAY   2164      //一日あたりの理想的な接種カロリー

struct Individual {
    double fitness;          // 適応度
    int chrom[CHROM_LENGTH]; // 染色体
    int order[CHROM_LENGTH];
};

typedef struct {
    int id;
    //char name[1024]
    double calorie;
    double protein;
    double fat;
    double carbohydrate;
}Product;

bool is_first_generation = true;

void load_products_data(char* fileName, Product* products) {

    int i, columnNum;
    char line[1024];
    FILE* fp;
    char* pos1;
    char* pos2;

    int dataNum;                    // データ数
    int colNum;                     // 項目数

    // 説明変数の数とデータ数の取得
    if ((fp = fopen(fileName, "r")) == NULL) {
        printf("%sが開けません\n", fileName);
        exit(1);
    }
    colNum = -1;
    dataNum = 0;
    while (fgets(line, 1024, fp)) {
        if (strcmp(line, "\n")) {
            columnNum = 1;
            pos1 = line;
            do {
                pos2 = strchr(pos1, ',');
                if (pos2) {
                    if (pos2 == pos1) {
                        printf("空のデータが含まれています．");
                        exit(1);
                    }
                    columnNum++;
                    pos1 = pos2 + 1;
                }
            } while (pos2);
            if (*pos1 == '\n') {
                printf("空のデータが含まれています．");
                exit(1);
            }
            if (colNum == -1) {
                colNum = columnNum - 1;
            }
            else if (colNum != columnNum - 1) {
                printf("列数の異なるレコードがあります．");
                exit(1);
            }
            dataNum++;
        }
    }
    fclose(fp);

    // データを読込む
    if ((fp = fopen(fileName, "r")) == NULL) {
        printf("%sが開けません\n", fileName);
        exit(1);
    }
    
   // ヘッダーを読み飛ばす
    fgets(line, 1024, fp);
   
   for(i = 0; i < dataNum - 1; i++) {
        fgets(line, 1024, fp);
        pos1 = line;
        //printf("%s", line);
        // id
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        //printf("%s\n", pos1);
        products[i].id = atoi(pos1);
        pos1 = pos2 + 1;
      
        //calorie
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        products[i].calorie = atof(pos1);
        pos1 = pos2 + 1;

        //protein
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        products[i].protein = atof(pos1);
        pos1 = pos2 + 1;

        //fat
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        products[i].fat = atof(pos1);
        pos1 = pos2 + 1;

        //carbohydrate
        pos2 = strchr(pos1, '\n');
        *pos2 = '\0';
        products[i].carbohydrate = atof(pos1);
        pos1 = pos2 + 1;
    }
    fclose(fp);
}

// 順序表現で初期化
void init_chrom_order_rep(int* chrom) {
    for (int i = 0; i < CHROM_LENGTH; ++i) {
        chrom[i] = 1 + (rand() % (CHROM_LENGTH - i));
    }
}

// 順序表現から巡回順生成
void chrom2order(Individual& ind) {
    std::vector<int> order_list(CHROM_LENGTH), order(CHROM_LENGTH);
    std::iota(order_list.begin(), order_list.end(), 1);
    for (int i = 0; i < CHROM_LENGTH; i++) {
        int chrom_i = ind.chrom[i];
        order[i] = order_list[chrom_i - 1];  
        order_list.erase(order_list.begin() + chrom_i - 1);
    }
    std::copy(order.begin(), order.end(), ind.order);  
}



void one_point_crossover(int* p1, int*p2, int* c1, int* c2) {
    int crossover_point = rand() % CHROM_LENGTH;
    
    // 子1と子2の前半部分をコピー
    std::copy(p1, p1 + crossover_point, c1);
    std::copy(p2, p2 + crossover_point, c2);

    // 子1と子2の後半部分を交換
    std::copy(p2 + crossover_point, p2 + CHROM_LENGTH, c1 + crossover_point);
    std::copy(p1 + crossover_point, p1 + CHROM_LENGTH, c2 + crossover_point);

    //     std::cout << std::endl;
    // for (int i = 0; i < CHROM_LENGTH; i++){
    //     std::cout << p1[i] << ", ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < CHROM_LENGTH; i++){
    //     std::cout << p2[i] << ", ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < CHROM_LENGTH; i++){
    //     std::cout << c1[i] << ", ";
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < CHROM_LENGTH; i++){
    //     std::cout << c2[i] << ", ";
    // }
}

void two_point_crossover(int* p1, int* p2, int* c1, int* c2) {
    int point1 = rand() % CHROM_LENGTH;
    int point2 = rand() % CHROM_LENGTH;

    // point1とpoint2の順番を正しいものにする（point1 < point2にする）
    if (point1 > point2) {
        std::swap(point1, point2);
    }

    // point1までの部分をコピー
    std::copy(p1, p1 + point1, c1);
    std::copy(p2, p2 + point1, c2);

    // point1からpoint2までの部分を交換してコピー
    std::copy(p2 + point1, p2 + point2, c1 + point1);
    std::copy(p1 + point1, p1 + point2, c2 + point1);

    // point2以降の部分をコピー
    std::copy(p1 + point2, p1 + CHROM_LENGTH, c1 + point2);
    std::copy(p2 + point2, p2 + CHROM_LENGTH, c2 + point2);
}

void uniform_crossover(int* p1, int* p2, int* c1, int* c2) {
    for (int i = 0; i < CHROM_LENGTH; i++) {
        if ((double)rand() / RAND_MAX < 0.5) {
            c1[i] = p1[i];
            c2[i] = p2[i];
        } else {
            c1[i] = p2[i];
            c2[i] = p1[i];
        }
    }
}

// グループの切れ目を取得する関数
std::vector<int> get_group_breaks(Individual &ind, Product* products) {
    std::vector<int> breaks;
    double current_calorie = 0.0;
    for (int j = 0; j < CHROM_LENGTH; j++) {
        current_calorie += products[ind.order[j]].calorie;
        if (current_calorie >= CAL_LIMIT_PER_DAY) {
            breaks.push_back(j);
            current_calorie = 0.0;
        }
    }
    return breaks;
}

void group_based_crossover(Individual &parent1, Individual &parent2, Individual &child1, Individual &child2, Product* products) {
    std::vector<int> breaks_p1 = get_group_breaks(parent1, products);
    std::vector<int> breaks_p2 = get_group_breaks(parent2, products);

    int crossover_point = -1;
    for (int i = 0; i < breaks_p1.size() && crossover_point == -1; i++) {
        for (int j = 0; j < breaks_p2.size(); j++) {
            if (breaks_p1[i] == breaks_p2[j]) {
                crossover_point = breaks_p1[i];
                break;
            }
        }
    }

    if (crossover_point != -1) {
        std::copy(parent1.chrom, parent1.chrom + crossover_point, child1.chrom);
        std::copy(parent2.chrom, parent2.chrom + crossover_point, child2.chrom);
        std::copy(parent2.chrom + crossover_point, parent2.chrom + CHROM_LENGTH, child1.chrom + crossover_point);
        std::copy(parent1.chrom + crossover_point, parent1.chrom + CHROM_LENGTH, child2.chrom + crossover_point);
    } else {
        // グループの切れ目が一致しない場合、通常の1点交叉を使用する（または他の方法を選択）
        std::copy(parent1.chrom, parent1.chrom + CHROM_LENGTH, child1.chrom);
        std::copy(parent2.chrom, parent2.chrom + CHROM_LENGTH, child2.chrom);
        // one_point_crossover(parent1.chrom, parent2.chrom, child1.chrom, child2.chrom);
    }
}


double calc_total_calorie(Product* products) {
    double total_calorie;
    total_calorie = 0;
    for (int i = 0; i < CHROM_LENGTH; i++) {
        total_calorie += products[i].calorie;
    }
    return total_calorie;
}
int binary_to_gray(int num) {
    return num ^ (num >> 1);
}

int gray_to_binary(int gray) {
    int num = 0;
    for (; gray; gray >>= 1) {
        num ^= gray;
    }
    return num;
}

void mutate_gray_code(int* num, int max_value) {
    int bit_length = (int)log2(max_value) + 1;
    int bit_pos = rand() % bit_length;
    int gray = binary_to_gray(*num);
    gray ^= (1 << bit_pos);  // Flip the bit at position bit_pos
    *num = gray_to_binary(gray);
}

void mutate_for_oreder_rep(int* chrom) {

    // int i = rand() % CHROM_LENGTH;
    // int value = chrom[i];
    // mutate_gray_code(&value, CHROM_LENGTH - i);
    // chrom[i] = 1 + (value % (CHROM_LENGTH - i));  // Ensure the mutated value is within the required range

    for (int i = 0; i < CHROM_LENGTH; i++) {
        if ( (double)rand()/RAND_MAX < PROBABILITY ) {
            chrom[i] =  1 + (rand() % (CHROM_LENGTH - i));
        }
    }
    
    
}



// ランキング選択
int rankingSelect()
{
   int num, denom, r;

   denom=POP_SIZE*(POP_SIZE+1)/2;
   r = ((rand() << 16) + (rand() << 1) + (rand() % 2)) % denom + 1;
   for (num = POP_SIZE; 0 < num; num--){
      if(r<=num){
         break;
      }
      r-=num;
   }
   return POP_SIZE-num;
}

void selectTwoDifferentParents(Individual* population, int &parent1, int &parent2) {
    // 最初の親を選択
    parent1 = rankingSelect();

    do {
        // 2番目の親を選択
        parent2 = rankingSelect();
    } while (parent1 == parent2);  // 最初の親と2番目の親が同じ場合は、再度選択
}

int rouletteSelect(Individual* population) {
    double total_fitness = 0.0;
    for (int i = 0; i < POP_SIZE; i++) {
        total_fitness += population[i].fitness;
    }

    double random_value = (double)rand() / RAND_MAX * total_fitness;
    double cumulative_fitness = 0.0;

    for (int i = 0; i < POP_SIZE; i++) {
        cumulative_fitness += population[i].fitness;
        if (cumulative_fitness >= random_value) {
            return i;
        }
    }

    // 通常はここに到達しないが、念のためのコード
    return POP_SIZE - 1;
}



// トーナメント選択
int tournamentSelect(Individual* population, int tournamentSize, int exclude = -1) {
    int best_index = -1;
    double best_fitness = -DBL_MAX;

    for (int i = 0; i < tournamentSize; i++) {
        int idx = rand() % POP_SIZE;

        // 排除するインデックスが指定されている場合
        while (idx == exclude) {
            idx = rand() % POP_SIZE;
        }

        if (population[idx].fitness > best_fitness) {
            best_fitness = population[idx].fitness;
            best_index = idx;
        }
    }

    return best_index;
}

void new_generation(Individual* population,  Product* products) {
    int i;
    int parent_index1, parent_index2;

    // mallocでIndividualの配列を動的に確保
    Individual* next_pop = (Individual*)malloc(POP_SIZE * sizeof(Individual));

    for(i = 0; i < POP_SIZE; i++) {
        next_pop[i].fitness = 0.0;
        memset(next_pop[i].chrom, 0, sizeof(next_pop[i].chrom));
        memset(next_pop[i].order, 0, sizeof(next_pop[i].order));
    }

    parent_index1 = 0;
    parent_index2 = 1;

    if(is_first_generation) {
        is_first_generation = false;
        init_chrom_order_rep(next_pop[0].chrom);
    } else {
        std::copy(population[0].chrom, population[0].chrom + CHROM_LENGTH, next_pop[0].chrom);
    }

    for (i = ELITE; i < POP_SIZE-1; i += 2) {
        //親選択
        selectTwoDifferentParents(population, parent_index1, parent_index2);
        //交叉
        group_based_crossover(population[parent_index1], population[parent_index2], next_pop[i], next_pop[i + 1], products);
        // one_point_crossover(population[parent_index1].chrom, population[parent_index2].chrom, next_pop[i].chrom, next_pop[i + 1].chrom);
    }   
    // std::cout << "nextpop" << std::endl;
    // for (int j = 0; j < CHROM_LENGTH; j++) {
    //     std::cout << next_pop[0].chrom[j] << ", ";
    // }
    

    //交配によって作る個体の数が奇数の時、最後の余った枠はランダムに生成した個体を入れる
    if (i != POP_SIZE) {
        init_chrom_order_rep(next_pop[i].chrom);
    }

    //突然変異
    for (i = ELITE; i < POP_SIZE; i++) {
            mutate_for_oreder_rep(next_pop[i].chrom);
    }

    for (int i = 0; i < POP_SIZE; i++) {
        chrom2order(next_pop[i]);
    }

    // コピー
    for(i = 0; i < POP_SIZE; i++) {

        population[i] = next_pop[i];

    }

    // 動的確保した配列を解放
    free(next_pop);
}




// 配列population[lb]～population[ub]をfitnessの昇順に整列（クイックソート）
// lb: 整列する範囲の下限の添字
// ub: 整列する範囲の上限の添字
void sort(Individual* population, int lb, int ub)
{
	int i, j, k;
	double pivot;
	Individual tmp;

	if(lb < ub) {
		k = (lb + ub) / 2;
		pivot = population[k].fitness;
		i = lb;
		j = ub;
		do {
			while(population[i].fitness < pivot)
				i++;
			while(population[j].fitness > pivot)
				j--;
			if(i <= j) {
				tmp = population[i];
				population[i] = population[j];
				population[j] = tmp;
				i++;
				j--;
			}
		} while(i <= j);
		sort(population, lb, j);
		sort(population, i, ub);
	}
}

//評価関数
void evaluation(Product* products, Individual* population) {
    int i, j;
    int term; //日数
    int product_index;
    //タンパク質、脂質、炭水化物、エネルギーの各日における合計量
    double day_protein, day_fat, day_carb, day_calorie;
    //タンパク質、脂質、炭水化物、エネルギー、各日の逸脱度
    double protein_deviance, fat_deviance, carb_deviance, calorie_deviance, day_deviance;
    double ave, sd, sqrSum;
    for (i = 0; i < POP_SIZE; i++) {
        j = 0;
        term = 0;
        day_protein = day_fat = day_carb = day_calorie = 0.0;
        ave = 0.0;
        sqrSum = 0.0;
        while (j < CHROM_LENGTH) {
            product_index = population[i].order[j];
            // printf("id:%d\n", products[product_index].id);
            // printf("protein:%f\n", products[product_index].protein);

            if (day_calorie < CAL_LIMIT_PER_DAY) {
                day_protein += products[product_index].protein;
                day_fat += products[product_index].fat;
                day_carb += products[product_index].carbohydrate;
                day_calorie += products[product_index].calorie;
                j++;
            }
            else {
                protein_deviance = fabs(0.2 - (double)(day_protein * 4) / day_calorie) * (1 / 0.2);
                fat_deviance = fabs(0.2 - (double)(day_fat * 9) / day_calorie) * (1 / 0.2);
                carb_deviance = fabs(0.6 - (double)(day_carb * 4) / day_calorie) * (1 / 0.6);
                calorie_deviance = fabs(IDEAL_CAL_PER_DAY - day_calorie) * (1 / IDEAL_CAL_PER_DAY);
                day_deviance = exp(protein_deviance) + exp(fat_deviance) + exp(carb_deviance) + exp(calorie_deviance) - 4;

                //最後に平均値と標準偏差を求めるために足しておく
                ave += day_deviance;
                sqrSum += day_deviance * day_deviance;
                term++;
                // printf("ave:%f\n", ave);
                // printf("sqrSum:%f\n", sqrSum);

                day_calorie = 0.0;
                day_protein = 0.0;
                day_fat = 0.0;
                day_carb = 0.0;
            }
        }
        //最終日の逸脱度を計算
        if (day_calorie != 0.0) {
            protein_deviance = fabs(0.2 - (double)(day_protein * 4) / day_calorie) * (1 / 0.2);
            fat_deviance = fabs(0.2 - (double)(day_fat * 9) / day_calorie) * (1 / 0.2);
            carb_deviance = fabs(0.6 - (double)(day_carb * 4) / day_calorie) * (1 / 0.6);
            calorie_deviance = fabs(IDEAL_CAL_PER_DAY - day_calorie) * (1 / IDEAL_CAL_PER_DAY);
            day_deviance = exp(protein_deviance) + exp(fat_deviance) + exp(carb_deviance) + exp(calorie_deviance) - 4;

            ave += day_deviance;
            sqrSum += day_deviance * day_deviance;
            term++;
        }

        ave /= term;                           //逸脱度の平均値
        sd = sqrt(sqrSum / term - ave * ave);  //逸脱度の標準偏差
        //printf("term:%d\n", term);
        //printf("innerdef_ave:%f\n", ave);
        //printf("innerdef_sqrSum:%f\n", sqrSum);
        //printf("innerdef_sd:%f\n", sd);
        population[i].fitness = ave + sd;
        //printf("innerdef_fitness:%f\n", population[i].fitness);
    }
    sort(population, 0, POP_SIZE - 1);
}

void save_result(char *fileName, int* chrom){
	FILE *fp;
	int i;
	// 書込み
	if((fp = fopen(fileName, "w")) == NULL) {
		printf("ファイル%sが開けません\n", fileName);
		exit(-1);
	}
	for(i = 0; i < CHROM_LENGTH; i++) {
		fprintf(fp, "%d,", chrom[i]);
	}
   fprintf(fp, "\n");
	fclose(fp);
}

int main() {
    int i;
    int gen;
    Product products[CHROM_LENGTH];
	char fname[64];

    Individual* population; 
    population = (Individual *)malloc(POP_SIZE * sizeof(Individual));

    // 乱数のタネの設定
    srand((unsigned int)time(NULL));

    //データをproductsに読み込む
	strcpy(fname, "dev_selected_500.csv");
    load_products_data(fname, products);

    // 初期集団生成
    for (i = 0; i < POP_SIZE; i++) {
        init_chrom_order_rep(population[i].chrom);
        chrom2order(population[i]);
    }


    // 進化
    for (gen = 1; gen <= GENERATION_MAX; gen++) {
        // 次世代生成
        new_generation(population,  products);
        // 評価
        evaluation(products, population);

        

        // for (int i = 0; i < POP_SIZE; i++)
        // {
        //     std::cout << population[i].fitness << std::endl;
        //     // std::copy(population[i].chrom, population[i].chrom + CHROM_LENGTH, std::ostream_iterator<int>(std::cout, ", "));
        //     std::copy(population[i].order, population[i].order + CHROM_LENGTH, std::ostream_iterator<int>(std::cout, ", "));
        // }
        
        // 途中経過表示
        fprintf(stderr, "第%4d世代\t最優良個体の適応度は%.16lf\n", gen, population[0].fitness);
    }
	strcpy(fname, "result.csv");
    save_result(fname, population[0].order);
    free(population);
    return 0;
}