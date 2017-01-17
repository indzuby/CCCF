#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <stdlib.h>
#include <time.h>
#pragma warning(disable:4996)

class Rating{
private:
	int user;
	int item;
	double rating;
public:
	Rating(int _user, int _item, double _rating) {
		user = _user;
		item = _item;
		rating = _rating;
	}
	int getUser() {
		return user;
	}
	int getItem() {
		return item;
	}
	double getRating() {
		return rating;
	}
	void setRating(double _rating) {
		rating = _rating;
	}
};



namespace distribution {

	template <typename RealType = double>
	class beta
	{
	public:
		typedef RealType result_type;

		class param_type
		{
		public:
			typedef beta distribution_type;

			explicit param_type(RealType a = 2.0, RealType b = 2.0)
				: a_param(a), b_param(b) { }

			RealType a() const { return a_param; }
			RealType b() const { return b_param; }

			bool operator==(const param_type& other) const
			{
				return (a_param == other.a_param &&
					b_param == other.b_param);
			}

			bool operator!=(const param_type& other) const
			{
				return !(*this == other);
			}

		private:
			RealType a_param, b_param;
		};

		explicit beta(RealType a = 2.0, RealType b = 2.0)
			: a_gamma(a), b_gamma(b) { }
		explicit beta(const param_type& param)
			: a_gamma(param.a()), b_gamma(param.b()) { }

		void reset() { }

		param_type param() const
		{
			return param_type(a(), b());
		}

		void param(const param_type& param)
		{
			a_gamma = gamma_dist_type(param.a());
			b_gamma = gamma_dist_type(param.b());
		}

		template <typename URNG>
		result_type operator()(URNG& engine)
		{
			return generate(engine, a_gamma, b_gamma);
		}

		template <typename URNG>
		result_type operator()(URNG& engine, const param_type& param)
		{
			gamma_dist_type a_param_gamma(param.a()),
				b_param_gamma(param.b());
			return generate(engine, a_param_gamma, b_param_gamma);
		}

		result_type min() const { return 0.0; }
		result_type max() const { return 1.0; }

		result_type a() const { return a_gamma.alpha(); }
		result_type b() const { return b_gamma.alpha(); }

		bool operator==(const beta<result_type>& other) const
		{
			return (param() == other.param() &&
				a_gamma == other.a_gamma &&
				b_gamma == other.b_gamma);
		}

		bool operator!=(const beta<result_type>& other) const
		{
			return !(*this == other);
		}

	private:
		typedef std::gamma_distribution<result_type> gamma_dist_type;

		gamma_dist_type a_gamma, b_gamma;

		template <typename URNG>
		result_type generate(URNG& engine,
			gamma_dist_type& x_gamma,
			gamma_dist_type& y_gamma)
		{
			result_type x = x_gamma(engine);
			return x / (x + y_gamma(engine));
		}
	};

	template <typename CharT, typename RealType>
	std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os,
		const beta<RealType>& beta)
	{
		os << "~Beta(" << beta.a() << "," << beta.b() << ")";
		return os;
	}

	template <typename CharT, typename RealType>
	std::basic_istream<CharT>& operator >> (std::basic_istream<CharT>& is,
		beta<RealType>& beta)
	{
		std::string str;
		RealType a, b;
		if (std::getline(is, str, '(') && str == "~Beta" &&
			is >> a && is.get() == ',' && is >> b && is.get() == ')') {
			beta = beta<RealType>(a, b);
		}
		else {
			is.setstate(std::ios::failbit);
		}
		return is;
	}


}


using namespace std;
using namespace distribution;

#define USER 943
#define ITEM 1682
#define RATINGS 100000
#define K 10
#define alphak1  1
#define alphak2  1
#define beta1 10
#define beta2 1

vector<Rating> rating_list;
vector<Rating> missing_list;

double y[USER+1][ITEM+1];
double pi1[USER + 1][K + 1];
double pi2[ITEM + 1][K + 1];
double theta[K + 1];
double pui[USER + 1][ITEM + 1];
double puik[USER + 1][ITEM + 1][K + 1];

double zuik[USER + 1][ITEM + 1][K+1];
double ziuk[ITEM + 1][USER + 1][K + 1];


void generativeCCCF(vector<Rating> O, vector<Rating> O_prime) {
	random_device rd;
	mt19937 gen(rd());
	beta<> beta(alphak1, alphak2), betak(beta1,beta2);
	bernoulli_distribution d;
	for (int i = 1; i <= K; i++) {
		for (int j = 1; j <= USER; j++)
			pi1[j][i] = beta(gen);
		for (int j = 1; j <= ITEM; j++)
			pi2[j][i] = beta(gen);
		theta[i] = betak(gen);
	}
	vector<Rating> merge = O;
	merge.insert(merge.end(), O_prime.begin(), O_prime.end());
	for (Rating r : merge) {
		double sum = 0;
		for (int i = 1; i <= K; i++) {
			d = bernoulli_distribution(pi1[r.getUser()][i]);
			zuik[r.getUser()][r.getItem()][i] = d(gen);

			d = bernoulli_distribution(pi2[r.getItem()][i]);
			zuik[r.getItem()][r.getUser()][i] = d(gen);


			if (zuik[r.getUser()][r.getItem()][i] == 1 && zuik[r.getItem()][r.getUser()][i] == 1) {
				puik[r.getUser()][r.getItem()][i] = theta[i];
				sum += (1 - puik[r.getUser()][r.getItem()][i]);
			}
			else
				puik[r.getUser()][r.getItem()][i] = 0;
		}
		pui[r.getUser()][r.getItem()] = 1 - sum;
		d = bernoulli_distribution(pui[r.getUser()][r.getItem()]);
		y[r.getUser()][r.getItem()] = d(gen);
		r.setRating(y[r.getUser()][r.getItem()]);
	}
}
void inferenceCCCF(vector<Rating> O, vector<Rating> O_prime) {
		
	int index = 1;
	srand(time(NULL));
	for (int i = 1; i <= K; i++)
		theta[i] = 1;
	for (int t = 0; t < 100; t++) {
		vector<Rating> merge = O;
		int size = O_prime.size();
		for (int i = 0;i<100 &&  index < size; index++,i++) {
			merge.push_back(O_prime[index]);
		}
		for (Rating rating : merge) {
			for (int k = 1; k <= K; k++) {
					// equation 7
			}
		}
		for (int tau = 0; tau < 100; tau++) { //converged check 
			for (int k = 1; k <= K; k++) {
				//equation 8
			}
		}
	}
	for (int i = 1; i <= USER; i++) {
		double sum = 0;
		int cnt = 0;
		for (int k = 1; k <= K; k++) {
			for (int j = 1; j <= ITEM; j++) {
				sum += zuik[i][j][k];
				if (zuik[i][j][k] > 0) cnt++;
			}
			pi1[i][k] = (sum + alphak1) / (cnt + alphak1 + alphak2);		//equation 9
		}
	}
	for (int i = 1; i <= ITEM; i++) {
		double sum = 0;
		int cnt = 0;
		for (int k = 1; k <= K; k++) {
			for (int j = 1; j <= USER; j++) {
				sum += ziuk[i][j][k];
				if (ziuk[i][j][k] > 0) cnt++;
			}
			pi2[i][k] = (sum + alphak1) / (cnt + alphak1 + alphak2);		//equation 10
		}
	}
}


void input() {
	FILE * in = fopen("u.data", "r");

	while (!feof(in)) {
		int user, item,rating,timestamp;
		fscanf(in, "%d %d %d %d", &user, &item, &rating, &timestamp);
		rating_list.push_back(Rating(user, item, rating/5.0));				// scale rating 0 ~ 1
		y[user][item] = rating / 5.0;
	}

	for (int i = 1; i <= USER; i++) {
		for(int j = 1; j<=ITEM;j++)
			if (y[i][j] <= 0) {
				missing_list.push_back(Rating(i, j, 0.0));
			}
	}

	fclose(in);
}

int main() {
	input();
	generativeCCCF(rating_list,missing_list);
	return 0;
}