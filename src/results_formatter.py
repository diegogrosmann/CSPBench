import statistics
from typing import List, Dict, Any
from tabulate import tabulate
from pathlib import Path

class ResultsFormatter:
    def __init__(self):
        self.results = {}
    
    def add_algorithm_results(self, algorithm_name: str, executions: List[Dict[str, Any]]):
        """Adiciona resultados de um algoritmo"""
        self.results[algorithm_name] = executions
    
    def format_detailed_results(self) -> str:
        """Formata todos os resultados com tabelas detalhadas"""
        output = []
        output.append("=" * 80)
        output.append("RELATÓRIO DETALHADO DE RESULTADOS")
        output.append("=" * 80)
        
        # Tabelas individuais por algoritmo
        for algorithm_name in self.results:
            output.append(self._format_algorithm_table(algorithm_name))
            output.append(self._format_algorithm_statistics(algorithm_name))
            output.append(self._format_algorithm_strings(algorithm_name))
            output.append("\n" + "=" * 80 + "\n")
        
        # Tabela comparativa final
        output.append(self._format_comparative_table())
        
        return "\n".join(output)
    
    def _format_algorithm_table(self, algorithm_name: str) -> str:
        """Formata tabela individual de um algoritmo"""
        executions = self.results[algorithm_name]
        
        output = [f"\n📊 TABELA DE EXECUÇÕES - {algorithm_name.upper()}"]
        output.append("-" * 60)
        
        # Cabeçalhos da tabela
        headers = [
            "Execução",
            "Tempo (s)",
            "Iterações",
            "Distância",
            "Solução Ótima",
            "String Encontrada"
        ]
        
        # Dados da tabela
        table_data = []
        for i, exec_data in enumerate(executions, 1):
            distancia = exec_data.get('distancia', exec_data.get('melhor_distancia', '-'))
            otima = "✅ Sim" if distancia == 0 else "❌ Não"
            string_result = exec_data.get('melhor_string', exec_data.get('string_resultado', ''))
            
            row = [
                i,
                f"{exec_data['tempo']:.4f}",
                exec_data.get('iteracoes', exec_data.get('num_iteracoes', 0)),
                distancia,
                otima,
                string_result[:20] + "..." if len(string_result) > 20 else string_result
            ]
            table_data.append(row)
        
        table = tabulate(table_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        
        return "\n".join(output)
    
    def _format_algorithm_statistics(self, algorithm_name: str) -> str:
        """Formata estatísticas detalhadas de um algoritmo"""
        executions = self.results[algorithm_name]
        
        # Extrair dados para análise com interface padronizada
        tempos = [exec_data['tempo'] for exec_data in executions]
        iteracoes = [exec_data.get('iteracoes', exec_data.get('num_iteracoes', 0)) for exec_data in executions]
        distancias = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) 
                     for exec_data in executions if exec_data.get('distancia', exec_data.get('melhor_distancia')) != '-']
        solucoes_otimas = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) == 0 
                          for exec_data in executions if exec_data.get('distancia', exec_data.get('melhor_distancia')) != '-']
        
        output = [f"\n📈 ESTATÍSTICAS DETALHADAS - {algorithm_name.upper()}"]
        output.append("-" * 60)
        
        # Estatísticas de tempo
        stats_data = [
            ["TEMPO DE EXECUÇÃO", ""],
            ["Média", f"{statistics.mean(tempos):.4f} s"],
            ["Mediana", f"{statistics.median(tempos):.4f} s"],
            ["Desvio Padrão", f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f} s"],
            ["Mínimo", f"{min(tempos):.4f} s"],
            ["Máximo", f"{max(tempos):.4f} s"],
            ["", ""],
            ["ITERAÇÕES", ""],
            ["Média", f"{statistics.mean(iteracoes):.2f}"],
            ["Mediana", f"{statistics.median(iteracoes):.2f}"],
            ["Desvio Padrão", f"{statistics.stdev(iteracoes) if len(iteracoes) > 1 else 0:.2f}"],
            ["Mínimo", f"{min(iteracoes)}"],
            ["Máximo", f"{max(iteracoes)}"],
        ]
        
        if distancias:
            stats_data.extend([
                ["", ""],
                ["DISTÂNCIA (QUALIDADE)", ""],
                ["Média", f"{statistics.mean(distancias):.2f}"],
                ["Mediana", f"{statistics.median(distancias):.2f}"],
                ["Desvio Padrão", f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}"],
                ["Melhor (Mínima)", f"{min(distancias)}"],
                ["Pior (Máxima)", f"{max(distancias)}"],
                ["", ""],
                ["DESEMPENHO", ""],
                ["Taxa de Sucesso", f"{sum(solucoes_otimas)}/{len(solucoes_otimas)} ({sum(solucoes_otimas)/len(solucoes_otimas)*100:.1f}%)"],
                ["Soluções Ótimas", f"{sum(solucoes_otimas)}"],
                ["Soluções Sub-ótimas", f"{len(solucoes_otimas) - sum(solucoes_otimas)}"],
            ])
        
        table = tabulate(stats_data, headers=["Métrica", "Valor"], tablefmt="grid", stralign="left")
        output.append(table)
        
        return "\n".join(output)
    
    def _format_algorithm_strings(self, algorithm_name: str) -> str:
        """Formata strings encontradas para auditoria"""
        executions = self.results[algorithm_name]
        
        output = [f"\n🔍 STRINGS PARA AUDITORIA - {algorithm_name.upper()}"]
        output.append("-" * 60)
        
        for i, exec_data in enumerate(executions, 1):
            distancia = exec_data.get('distancia', exec_data.get('melhor_distancia', '-'))
            string_result = exec_data.get('melhor_string', exec_data.get('string_resultado', ''))
            iteracoes = exec_data.get('iteracoes', exec_data.get('num_iteracoes', 0))
            
            status = "✅ SOLUÇÃO ÓTIMA" if distancia == 0 else "⚠️ SOLUÇÃO SUB-ÓTIMA"
            output.append(f"Execução {i:2d} ({status}):")
            output.append(f"  String: '{string_result}'")
            output.append(f"  Distância: {distancia}")
            output.append(f"  Iterações: {iteracoes}")
            output.append(f"  Tempo: {exec_data['tempo']:.4f}s")
            output.append("")
        
        return "\n".join(output)
    
    def _format_comparative_table(self) -> str:
        """Formata tabela comparativa entre algoritmos"""
        output = ["\n🏆 TABELA COMPARATIVA FINAL"]
        output.append("=" * 80)
        
        if not self.results:
            output.append("Nenhum resultado disponível para comparação.")
            return "\n".join(output)
        
        # Cabeçalhos
        headers = [
            "Algoritmo",
            "Tempo Médio (s)",
            "Iterações Médias",
            "Melhor Distância",
            "Distância Média",
            "Taxa Sucesso (%)",
            "Desvio Tempo"
        ]
        
        # Calcular métricas comparativas
        comparative_data = []
        for algorithm_name, executions in self.results.items():
            tempos = [exec_data['tempo'] for exec_data in executions]
            iteracoes = [exec_data.get('iteracoes', exec_data.get('num_iteracoes', 0)) for exec_data in executions]
            distancias = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) 
                         for exec_data in executions if exec_data.get('distancia', exec_data.get('melhor_distancia')) != '-']
            solucoes_otimas = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) == 0 
                              for exec_data in executions if exec_data.get('distancia', exec_data.get('melhor_distancia')) != '-']
            
            if distancias:
                row = [
                    algorithm_name,
                    f"{statistics.mean(tempos):.4f}",
                    f"{statistics.mean(iteracoes):.1f}",
                    f"{min(distancias)}",
                    f"{statistics.mean(distancias):.2f}",
                    f"{sum(solucoes_otimas)/len(solucoes_otimas)*100:.1f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}"
                ]
            else:
                row = [
                    algorithm_name,
                    f"{statistics.mean(tempos):.4f}",
                    f"{statistics.mean(iteracoes):.1f}",
                    "-",
                    "-",
                    "0.0",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}"
                ]
            comparative_data.append(row)
        
        # Ordenar por melhor distância (crescente), depois por tempo médio
        comparative_data.sort(key=lambda x: (float('inf') if x[3] == '-' else float(x[3]), float(x[1])))
        
        table = tabulate(comparative_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        
        # Adicionar ranking
        output.append("\n🥇 RANKING POR PERFORMANCE:")
        output.append("-" * 40)
        for i, row in enumerate(comparative_data, 1):
            medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else f"{i}°"
            output.append(f"{medal} {row[0]} - Distância: {row[3]} | Tempo: {row[1]}s | Sucesso: {row[5]}%")
        
        return "\n".join(output)
    
    def save_detailed_report(self, filename: str = "relatorio_detalhado.txt"):
        """Salva relatório detalhado em arquivo"""
        results_dir = Path(__file__).parent.parent / "results"
        results_dir.mkdir(exist_ok=True)
        filepath = results_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(self.format_detailed_results())
