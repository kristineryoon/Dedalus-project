# DEDALUS 튜토리얼 2

## OVERVIEW:

-Dedalus에서 Field와 연산자 객체의 상호작용과 기본적인 설정 방법을 익힌다. 

-Dedalus는 수학적 표현식과 PDE을 표현하기 위해 필드와 연산자 추상화를 사용하여 기호 대수 시스템을 설정한다.

## **2.1: Fields**

**Creating a field**

-Field 객체는 기저 집합 위에 정의된 스칼라값  필드를 나타낸다.

-여기서 Field는 분배기, 기저들의 리스트, 선택적으로 이름을 전달하여 직접 Field class로 정의 가능하다.

```
coords = d3.CartesianCoordinates('x', 'y')
dist = d3.Distributor(coords, dtype=np.float64)

xbasis = d3.RealFourier(coords['x'], 64, bounds=(-np.pi, np.pi), dealias=3/2)
ybasis = d3.Chebyshev(coords['y'], 64, bounds=(-1, 1), dealias=3/2)

f = dist.Field(name='f', bases=(xbasis, ybasis))
```

**Vector and tensor fields**

-Field 클래스는 스칼라 값 필드를 생성

-벡터값 필드는 Vector Field를 통해 만들어지는데 벡터가 어떤 성분 방향을 가질지에 해당하는 좌표계를 함께 전달해야한다.

-이것은 해당좌표계의 접다발을 필드 다발로 지정을 의미

-임의 차수 텐서 필드도 Tensor Field  생성자를 사용하고 좌표계들은 튜플을 전달 함으로서 생성 할수 있다.(텐서 다발 기술)

주의점

필드의 기저는 그 필드의 공간적 변화를 나타내는 반면, 벡터나 텐서 다발은 필드의 성분을 나타낸다는 점을 기억하자.

예를 들어, x 성분과 y 성분을 가진 2차원 벡터가 있다고 하면 이벡터가 오직 x 방향으로만 변하고 y 방향은 변하지 않는다면 이 벡터 필드는 x 방향에 대한 기저만 갖는다.

**Manipulating field data**

Field 객체에는 데이터를 서로다른  Layout의 형태에 따라 변환할수 있는 다양한 method들이 존재 

필드 데이터는 어떤 레이아웃이든, 해당 레이아웃 객체를 이용하여 필드를 인덱싱하거나 불러올수 있다.

병렬 환경에서 필드 데이터를 접근할 때, 각 프로세스는 전역적으로 분리된 데이터셋 중 자신에게 할당된 지역 데이터만 다룬다.

즉, 도메인 객체가 제공하는 지역 그리드를 사용하면, 필드 그리드 데이터를 병렬 안전한 방식으로 설정 가능

```python
x = dist.local_grid(xbasis)
y = dist.local_grid(ybasis)
f['g'] = np.exp((1-y**2)*np.cos(x+np.cos(x)*y**2)) * (1 + 0.05*np.cos(10*(x+2*y)))

# Plot grid values
plot_bot_2d(f, figkw=figkw, title="f['g']");
```

![image.png](image.png)

스펙트럴 계수공간으로 변환

```
f['c']

# Plot log magnitude of spectral coefficients
log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
plot_bot_2d(f, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])", figkw=figkw)
```

![image.png](image%201.png)

필드의 스펙트럴 계수를 살펴보는 것은 매우 유용하다. 가장 높은 모드의 진폭이 필드의 스펙트럴 이산화 과정에서 발생하는 절단 오차를 잘 보여준다.

**Vector and tensor components**

벡터장과 텐서장은 첫 번째 축에 그 성분들을 포함하는 더 높은 차원의 데이터 배열을 갖는다.

`u['g'].shape`

`(2, 64, 1)`

첫 번째 축의 길이는 2로, 이는 벡터의 두 성분에 해당한다.

이 벡터장은 z 방향으로는 변화가 없기 때문

field.shape = (components, Ny, Nx) 라는 관례

- 카르테시안 좌표계에서는 벡터/텐서 성분이 단순히 각 방향 성분으로 분리되어 해석이 쉽다.
- 하지만 구형 등 다른 좌표계에서는 변환 과정에서 성분이 섞이므로 해석이 복잡하다.
- 그래서 초기 데이터 설정은 **격자 공간**에서 하는 것이 좋다.

**Field scale factors**

필드 객체의 **change_scales** 메서드는 필드의 데이터를 격자 공간으로 변환할때 사용되는 스케일링 인자를 변경하는데 사용된다.

격자 배열을 사용해 필드 데이터를 설정 할때, 필드 스케일 인자와 격자 스케일인자가 일치하지 않으면 할상 shape 오류가 발생한다.

큰 스케일 인자-필드 데이터를 고해상도 격자로 보간 가능

작은 스케일 인자-필드의 저해상도 격자 표현을 볼수 있다.

### 🧠 요약 정리

| 개념 | 설명 |
| --- | --- |
| **`change_scales` 메서드** | 격자 공간 변환 시 스케일(해상도)을 조정하는 함수 |
| **큰 scale factor (>1)** | 고해상도 격자, 더 많은 점으로 보간 |
| **작은 scale factor (<1)** | 저해상도 격자, 데이터 손실 발생 가능 |
| **주의사항** | 필드와 격자의 스케일이 맞지 않으면 shape error 발생 |

```
f.change_scales(4)

# Plot grid values
f['g']
plot_bot_2d(f, title="f['g']", figkw=figkw);
```

![image.png](image%202.png)

## **2.2: Operators**

**Arithmetic with fields(필드간 사칙연산)**

필드에 대한 수학적 연산은 사칙연산 미분 적분 보간 등은 직접 수치계산을 수행하지 않고 operater class로 표현된다.

이 operator 객체는 연산을 수행하라라는 수학적 식 자체를 저장한다.

필드들끼리, 또는 필드와 스칼라(숫자) 사이의 산술 연산은 파이썬의 기본 산술 연산자(중위 연산자, infix operators) 를 그대로 사용해서 수행할 수 있습니다.

`g_op = 1 - 2*f
print(g_op)`

`C(C(1)) + -1*2*f`

이 계산 결과  얻는 객체는 다른 필드로 얻는것이 아니라 계산식에 해당하는 연산을 나타내는 연산자 객체

이 연산을 실제로 계산하려면 **evaluate()**라는 메서드를 사용해야 하며, 이 메서드는 계산 결과를 담은 새로운 필드를 반환한다.

모든 연산자들의 계산에는 기저가 생성될때 설정된 dealias 스케일 인자가 사용된다.

```
g = g_op.evaluate()

# Plot grid values
g['g']
plot_bot_2d(g, title="g['g']", figkw=figkw);
```

![image.png](image%203.png)

**Building expressions**

연산자 객체는 다른 연산자의 인자로 전달 될수 있으며, 여러 연산이 트리 형태로 결합되어 더 복잡한 수학적 표현이 가능하다.

```
h_op = 1 / np.cosh(g_op + 2.5)
print(h_op)

`Pow(cosh(C(C(1)) + -1*2*f + C(C(2.5))), -1)`
```

POW는 역수

dedalus.tools 모듈에 있는 도움 도구(helper) 를 사용하면 이 연산자의 구조(operator’s structure) 를 시각적으로(그래프로) 그릴 수 있습니다.

```
from dedalus.tools.plot_op import plot_operator
plot_operator(h_op, figsize=6, fontsize=14, opsize=0.4)
```

![image.png](image%204.png)

```
h = h_op.evaluate()

# Plot grid values
h['g']
plot_bot_2d(h, title="h['g']", figkw=figkw);
```

![image.png](image%205.png)

**Deferred evaluation(지연 평가)**

연산자 객체는 필드들에 대한 연산을 기호적으로 표현

이 연산자들은 지연 평가 방식을 사용해 계산

만약 연산에 사용된 필드의 데이터를 변경한 다음, 연산자를 다시 평가하면 새로운 결과를 얻을 수 있다.

```
# Change scales back to 1 to build new grid data
f.change_scales(1)
f['g'] = 3*np.cos(1.5*np.pi*y)**2 * np.cos(x/2)**4 + 3*np.exp(-((2*x+2)**2 + (4*y+4/3)**2)) + 3*np.exp(-((2*x+2)**2 + (4*y-4/3)**2))

# Plot grid values
f['g']
plot_bot_2d(f, title="f['g']", figkw=figkw);
```

![image.png](image%206.png)

```
h = h_op.evaluate()

# Plot grid values
h['g']
plot_bot_2d(h, title="h['g']", figkw=figkw);
```

![image.png](image%207.png)

즉, `f`를 재정의하면, 그것을 참조하는 모든 연산자(`g_op`, `h_op` 등)의 계산 결과가 전부 달라집니다.

연관 되는 데이터 다 달라짐!!

**Differential operators**

연산자는 미분 필드에서 사용된다. 부분 미분 연산자는 한 차원의 기저에 대해 미분 연산자를 사용하고, 미분할 좌표를 지정 함으로서 구
`fx = d3.Differentiate(f, coords['x'])`

다차원의 경우 

1차원은 내장된 미분 연산을 사용하지만, 2차원 이상에서는 

- `Gradient` for arbitrary fields.
- `Divergence` for arbitrary vector/tensor fields.
- `Curl` for vector fields.
- `Laplacian`, defined as the divergence of the gradient, for arbitrary fields.

**Tensor operators**

| 연산자 | 기능 | 기본 인덱스 | 설명 |
| --- | --- | --- | --- |
| **Trace** | 두 인덱스 수축 (대각합 계산) | 0, 1 | 텐서의 두 축을 합해 차원을 줄임 |
| **TransposeComponents** | 두 인덱스 교환 (전치) | 0, 1 | 지정된 두 축의 위치를 바꿈 |
| **Skew** | 2D 벡터를 90° 회전 | 해당 없음 | (x, y) → (−y, x) |

**Integrals and averages**

스칼라 필드(scalar field)의 **좌표 또는 좌표계에 대한 적분(integral)**과 **평균(average)**은 **Integrate** 및 **Average** 연산자를 사용하여 계산된다.

**Interpolation**

좌표를 따라 수행되는 **보간(interpolation)**은 **Interpolate 연산자(Interpolate operator)**를 사용하거나,

필드(field)나 연산자(operator)의 **`__call__` 메서드**를 이용하여 수행할 수 있다.

이때, **좌표(coordinate)**와 **위치(position)**는 키워드(keyword) 인자로 지정한다.

편의를 위해, 문자열 **'left'**, **'right'**, **'center'**를 사용하여

1차원 구간(1D interval)의 **왼쪽 끝점**, **오른쪽 끝점**, 그리고 **중간점**을 각각 가리킬 수 있다.
