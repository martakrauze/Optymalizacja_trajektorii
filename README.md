# Optymalizacja_trajektorii
Praca inżynierska - Optymalizajca trajektorii płaskiego układu wieloczłonowego metodą bezpośredniej transkrypcji.

Skrypty napisane dla GNU Octave.

## Przykładowe uruchomienie skryptu

```
$ octave three_parts_opt.m
$ octave --persist show.m
```

## Opis folderów i skryptów:
### Drobne_przykłady
Własna implementacja przykładów z publikacji M. Kelley https://epubs.siam.org/doi/10.1137/16M1062569

-> block_optim.m - przykład optymalizacji trajektorii przesuwania klocka

-> cart_pole.m - przykład optymalizacji trajektorii wózka z wahadłem podczas manewru "swing-up"

### Układ_dwuczłonowy

Skrypty dla układu dwuczłonowego z parami kinematycznymi obrotowymi

-> two_parts_opt.m - skrypt do optymalizacji trajektorii

-> show.m - skrypt do wyświetlania wyników optymalizacji trajektorii

-> two_parts_dynamics_function.m - funkcja obliczająca dynamikę układu we współrzędnych złączowych

### Układ_trójczłonowy

Skrypty dla układu trójczłonowego z parami kinematycznymi obrotowymi

#### Dynamika

-> three_parts_dyn.m - skrypt do wyznaczania ruchu układu

-> three_parts_dynamics_function - funkcja obliczająca dynamikę układu we współrzędnych złączowych

-> *.txt - pliki zawieracjące wyniki z analizy w programie Freedyn (http://www.freedyn.at/) do porównania

#### Optymalizacja

-> three_parts_opt.m - skrypt do optymalizacji trajektorii z punktu do punktu (najmniejsze momenty sił)

-> three_parts_tvp.m - skrypt do optymalizacji trajektorii tak by końcówka układu poruszała się jak najbliżej trajektorii referencyjnej po linii prostej z trapezowym profilem prędkości

-> show.m - skrypt do wyświetlania wyników optymalizacji trajektorii

-> three_parts_dynamics_function.m - funkcja obliczająca dynamikę układu we współrzędnych złączowych

-> function_tvp.m - funkcja obliczająca trajektorię z trapezowym profilem prędkości

-> three_parts_opt_analiza.m - skrypt wykonujący obliczenia optymalizacyjne dla zmieniającej się od 3 do 40 liczby kroków czasowych metody trapezów (Uwaga: Czas obliczeń tego skryptu to kilka godzin)

-> show_por.m - skrypt wyświetlający porównanie optymalnej trajektorii dla różnej liczby kroków czasowych