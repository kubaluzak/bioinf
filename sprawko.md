# Rekonstrukcja sekwencji DNA na podstawie sond binarnych — opis problemu i algorytmu

---

## Autorzy  
**Jakub Urbaniak**  
**Wiktor Wachowski**

---

## Spis treści

0. [Opis i formalizacja problemu](#0-opis-i-formalizacja-problemu)  
1. [Wprowadzenie i opis problemu](#1-wprowadzenie-i-opis-problemu)  
2. [Algorytm BFS — dokładny](#2-algorytm-bfs--dokładny)  
   2.1. [Reprezentacja stanu](#21-reprezentacja-stanu)  
   2.2. [Inicjalizacja BFS](#22-inicjalizacja-bfs)  
   2.3. [Iteracyjny proces BFS](#23-iteracyjny-proces-bfs)  
   2.4. [Zakończenie i wynik](#24-zakończenie-i-wynik)  
3. [Analiza złożoności czasowej](#3-analiza-złożoności-czasowej)  
   3.1. [Parametry](#31-parametry)  
   3.2. [Przestrzeń stanów](#32-przestrzeń-stanów)  
   3.3. [Koszt pojedynczego stanu](#33-koszt-pojedynczego-stanu)  
   3.4. [Całkowita złożoność](#34-całkowita-złożoność)  
4. [Zależność czasu działania algorytmu od długości sekwencji](#4-zależność-czasu-działania-algorytmu-od-długości-sekwencji)

---

## 0. Opis i formalizacja problemu

Celem problemu jest rekonstrukcja pełnej sekwencji DNA na podstawie ograniczonych i zakodowanych danych sondowych. Dane sondy są generowane w dwóch różnych kodowaniach binarnych, zwanych spektrami: WS (Weak/Strong) oraz RY (Purine/Pyrimidine).

### Dane wejściowe

- **Sekwencja startowa** \( S \) o długości \( m \) — znany fragment DNA, od którego rozpoczynamy rekonstrukcję,  
- **Sondy** w dwóch spektrach:  
  - Spektrum WS: sondy składające się z liter \(\{W, S\}\) reprezentujące nukleotydy według reguł \(A, T \to W\) i \(C, G \to S\),  
  - Spektrum RY: sondy składające się z liter \(\{R, Y\}\) reprezentujące nukleotydy według reguł \(A, G \to R\) i \(C, T \to Y\),  
- Długość docelowej sekwencji \( n \), do której chcemy zrekonstruować pełną sekwencję.

### Charakterystyka problemu

Każda sonda ma długość \( k \) i jest zdefiniowana jako prefix długości \( k-1 \) oraz nukleotyd \( N \) dopisywany na końcu. Jednak z powodu podwójnego kodowania jedna sonda w spektrum WS może odpowiadać wielu różnym k-merom DNA, podobnie jak sonda w spektrum RY. Dzięki analizie obu spektrów jednocześnie można ograniczyć liczbę możliwych rozszerzeń i zredukować niejednoznaczności.

### Formalizacja

Niech \( \Sigma = \{A, C, G, T\} \) będzie alfabetem nukleotydów. Definiujemy dwie funkcje kodujące:

\[
\text{encode}_{WS}: \Sigma \to \{W, S\} \quad\text{gdzie}\quad A, T \mapsto W; \quad C, G \mapsto S
\]

\[
\text{encode}_{RY}: \Sigma \to \{R, Y\} \quad\text{gdzie}\quad A, G \mapsto R; \quad C, T \mapsto Y
\]

Dla aktualnego suffixu rekonstrukcji o długości \( k-1 \) w obu kodowaniach szukamy sond, których prefixy pasują do tych suffixów, a nukleotyd \( N \) dopisujemy tylko wtedy, gdy jest zgodny w obu spektrach.

### Cel

Znaleźć ciąg \( S' \in \Sigma^n \), który:

- Rozszerza znany start \( S \),  
- Jest spójny z dostępnymi sondami w obu spektrach (dla każdego dopisanego nukleotydu istnieją sondy w spektrum WS i RY o odpowiednich prefixach i dopisanych nukleotydach),  
- Minimalizuje błędy wynikające z niedokładności sond.

---

## 1. Wprowadzenie i opis problemu

Rekonstrukcja sekwencji DNA to zadanie odtworzenia pełnej sekwencji nukleotydów na podstawie fragmentarycznych informacji, takich jak krótkie fragmenty DNA, zwane sondami. W klasycznych metodach każda sonda reprezentuje jednoznacznie określony fragment sekwencji, co pozwala na prostą, liniową rekonstrukcję. Jednak w nowoczesnych technikach, w szczególności wykorzystujących kodowanie binarne dwóch spektrów, problem staje się bardziej złożony.

W naszym podejściu korzystamy z dwuspektralnego kodowania sond, gdzie każdy nukleotyd jest jednocześnie reprezentowany w dwóch niezależnych kodowaniach: w spektrum WS (gdzie A i T kodowane są jako "W", C i G jako "S") oraz w spektrum RY (gdzie A i G są purynami kodowanymi jako "R", a C i T są pirymidynami kodowanymi jako "Y"). Takie podwójne kodowanie zwiększa ilość informacji i zmniejsza liczbę błędów, jednak sondy nie są już jednoznaczne — pojedyncza sonda może odpowiadać wielu możliwym fragmentom oryginalnej sekwencji.

Zadaniem algorytmu jest odtworzenie pełnej sekwencji DNA o zadanej długości, rozpoczynając od znanego fragmentu startowego. Algorytm musi uwzględnić informacje z obu spektrów, tak aby dopasować sondy z pierwszego i drugiego spektrum jednocześnie, eliminując możliwe błędy i niezgodności.

---


## 2. Algorytm BFS — dokładny

Do rozwiązania problemu wykorzystujemy algorytm **przeszukiwania wszerz (BFS)** po przestrzeni stanów reprezentującej możliwe rozszerzenia rekonstrukcji sekwencji.

### 2.1. Reprezentacja stanu

Każdy stan algorytmu reprezentowany jest przez:

- **Suffix sekwencji** z ostatnich \(k-1\) nukleotydów zakodowanych w spektrum WS,  
- **Suffix sekwencji** z ostatnich \(k-1\) nukleotydów zakodowanych w spektrum RY,  
- Aktualnie zrekonstruowaną sekwencję nukleotydów (ciąg liter A, C, T, G).

W kodzie, zamiast przechowywać całe sekwencje, operujemy na suffixach i rozszerzamy je o kolejne nukleotydy zgodne z sondami.

### 2.2. Inicjalizacja BFS

- Na początku kodujemy suffix startowej sekwencji w obu spektrach,  
- Tworzymy dwie mapy prefixów (po jednym na każde spektrum), które dla prefixu długości \(k-1\) podają zbiór możliwych nukleotydów do dopisania,  
- Inicjujemy kolejkę BFS z początkowym stanem: suffixy startowe i startowa sekwencja,  
- Tworzymy zbiór `visited` by pamiętać, które pary suffixów już zostały przetworzone — to zapobiega powtórkom.

### 2.3. Iteracyjny proces BFS

W każdej iteracji:

1. Pobieramy stan z kolejki — aktualne suffixy WS i RY oraz dotychczas zrekonstruowaną sekwencję.  
2. Sprawdzamy, czy długość sekwencji osiągnęła wymaganą wartość — jeśli tak, to zwracamy rozwiązanie.  
3. Pobieramy z map prefixów kandydujące nukleotydy, które można dołączyć do aktualnych suffixów WS i RY.  
4. Iterujemy przez nukleotydy \( \{A, C, T, G\} \) i sprawdzamy, czy nukleotyd pasuje jednocześnie do obu spektrów (czy występuje jako dozwolony w obu mapach dla odpowiednich prefixów).  
5. Dla każdego pasującego nukleotydu tworzymy nowy stan:  
   - Obliczamy nowe suffixy, przesuwając istniejący suffix o 1 pozycję i dopisując kodowaną literę (w spektrum WS i RY),  
   - Doklejamy nukleotyd do aktualnej sekwencji,  
   - Jeśli para suffixów nie była wcześniej odwiedzona, dodajemy nowy stan do kolejki i do zbioru `visited`.

### 2.4. Zakończenie i wynik

- Algorytm działa dopóki kolejka nie jest pusta lub nie znajdzie sekwencji o wymaganej długości,  
- W przypadku braku rozwiązania informuje o tym użytkownika.

---

## 3. Analiza złożoności czasowej

### 3.1. Parametry

| Symbol | Znaczenie                       |
|--------|--------------------------------|
| \(k\)  | Długość sondy                  |
| \(n\)  | Długość docelowej sekwencji   |
| \(S\)  | Liczba sond w każdym spektrum |

### 3.2. Przestrzeń stanów

Algorytm operuje na parach suffixów długości \(k-1\) z obu spektrów.

Liczba stanów do rozpatrzenia:

\[
O(S^2)
\]

### 3.3. Koszt pojedynczego stanu

- Pobranie kandydatów i sprawdzenie odwiedzenia to operacje w czasie stałym \(O(1)\).

### 3.4. Całkowita złożoność

\[
O(S^2)
\]

z założeniem, że eksplorowane są wszystkie możliwe pary prefixów. W praktyce jest zwykle mniejsza dzięki filtrowaniu i zatrzymaniu po znalezieniu rozwiązania.

---

## 4. Zależność czasu działania algorytmu od długości sekwencji
Parametry dla każdego przypadku: 

k = 10

sqpe = n/4

pose = n/2

| Długość sekwencji | Czas wykonania [s] | Przykładowa sekwencja (fragment)               |
|-------------------|--------------------|-----------------------------------------------|
| 20                | 0.000004           | ATGTTGTAAACCAAGATCTA                          |
| 30                | 0.000026           | GGGCCTTGTCTTTGCCTACTGAAGATATAT                |
| 50                | 0.000063           | GGGAGTACAGGCTTGAGCATATTAAAAAGACTGCTTTTTAATAT |
| 100               | 0.000196           | CCATCTTGCTCTCTCGATGCTGCCAGGCGTAGAAGTTGAGC... |
| 150               | 0.000279           | CCCGTCTCTACTAAAAATACAAAATCAGCCAAGCGTGGTGGC... |
| 250               | 0.000481           | ATCTGGAGGGAAAACTGTGTGGAGAGAACACTTGACAAGAAA... |
| 500               | 0.001098           | AAGTGGGAGAAAAAGCTGCTGCCCATCCAGCAATGGAGCTTC... |
| 1000              | 0.002375           | AATCAACTCAAAATTGATTGAACATGTAAATATAAGACCTGA... |

---
