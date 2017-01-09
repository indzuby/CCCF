#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <vector>

class Rating{
private:
	int user;
	int item;
	int rating;
	int timestamp;
public:
	Rating(int _user, int _item, int _rating ,int _timestamp) {
		user = _user;
		item = _item;
		rating = _rating;
		timestamp = _timestamp;
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

double y[USER+1][ITEM+1];
int pi1[USER + 1][K + 1];
int pi2[ITEM + 1][K + 1];
int theta[K + 1];
int pui[USER + 1][ITEM + 1];
int puik[USER + 1][ITEM + 1][K + 1];

int zuik[USER + 1][ITEM + 1][K+1];
int ziuk[ITEM + 1][USER + 1][K + 1];



int main() {
	random_device rd;
	mt19937 gen(rd());
	beta<> beta(2, 2);
	bernoulli_distribution d;
	for (int i = 0; i < 100; i++) {
		double p = beta(gen);
		cout <<"beta : "<< p << endl;
		d = bernoulli_distribution(p);
		cout << "bernoulli : " << d(gen) << endl;
	}

	return 0;
}