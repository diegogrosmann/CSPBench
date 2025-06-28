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
            "Distância",
            "Status"
        ]
        
        # Dados da tabela
        table_data = []
        for i, exec_data in enumerate(executions, 1):
            distancia = exec_data.get('distancia', exec_data.get('melhor_distancia', '-'))
            status = "✓ OK"
            
            # Verifica se houve erro
            if exec_data.get('erro'):
                status = f"✗ {exec_data['erro']}"
                distancia = "-"
            elif exec_data.get('timeout'):
                status = "⏰ Timeout"
            elif distancia == float('inf'):
                status = "∞ Sem solução"
                distancia = "∞"
            
            row = [
                i,
                f"{exec_data['tempo']:.4f}",
                distancia,
                status
            ]
            table_data.append(row)
        
        table = tabulate(table_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        
        return "\n".join(output)
    
    def _format_algorithm_statistics(self, algorithm_name: str) -> str:
        """Formata estatísticas detalhadas de um algoritmo"""
        executions = self.results[algorithm_name]
        
        # Filtrar apenas execuções válidas
        valid_executions = [
            exec_data for exec_data in executions 
            if not exec_data.get('erro') and 
               exec_data.get('distancia', exec_data.get('melhor_distancia')) not in ['-', float('inf')]
        ]
        
        output = [f"\n📈 ESTATÍSTICAS DETALHADAS - {algorithm_name.upper()}"]
        output.append("-" * 60)
        
        if not valid_executions:
            output.append("❌ Nenhuma execução válida para calcular estatísticas.")
            return "\n".join(output)
        
        # Extrair dados para análise
        tempos = [exec_data['tempo'] for exec_data in valid_executions]
        distancias = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) 
                     for exec_data in valid_executions]
        
        # Estatísticas de tempo
        stats_data = [
            ["TEMPO DE EXECUÇÃO", ""],
            ["Execuções Válidas", f"{len(valid_executions)}/{len(executions)}"],
            ["Média", f"{statistics.mean(tempos):.4f} s"],
            ["Mediana", f"{statistics.median(tempos):.4f} s"],
            ["Desvio Padrão", f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f} s"],
            ["Mínimo", f"{min(tempos):.4f} s"],
            ["Máximo", f"{max(tempos):.4f} s"],
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
            
            output.append(f"Execução {i:2d}:")
            if exec_data.get('erro'):
                output.append(f"  ❌ Erro: {exec_data['erro']}")
                output.append(f"  Tempo: {exec_data['tempo']:.4f}s")
            else:
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
            "Desvio Tempo",
            "Melhor Distância",
            "Distância Média", 
            "Desvio Distancia",
            "Taxa Sucesso (%)"
        ]
        
        # Calcular métricas comparativas
        comparative_data = []
        for algorithm_name, executions in self.results.items():
            # Filtrar execuções válidas
            valid_executions = [
                exec_data for exec_data in executions 
                if not exec_data.get('erro') and 
                   exec_data.get('distancia', exec_data.get('melhor_distancia')) not in ['-', float('inf')]
            ]
            
            if not valid_executions:
                # Algoritmo falhou em todas as execuções
                tempos = [exec_data['tempo'] for exec_data in executions]
                row = [
                    algorithm_name,
                    f"{statistics.mean(tempos):.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    "ERRO",
                    "ERRO", 
                    "ERRO",
                    "0.0"
                ]
            else:
                tempos = [exec_data['tempo'] for exec_data in valid_executions]
                distancias = [exec_data.get('distancia', exec_data.get('melhor_distancia', float('inf'))) 
                             for exec_data in valid_executions]
                
                # Taxa de sucesso: execuções válidas / total de execuções
                taxa_sucesso = len(valid_executions) / len(executions) * 100
                solucoes_otimas = sum(1 for d in distancias if d == 0)
                taxa_otima = solucoes_otimas / len(distancias) * 100 if distancias else 0
                
                row = [
                    algorithm_name,
                    f"{statistics.mean(tempos):.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    f"{min(distancias)}",
                    f"{statistics.mean(distancias):.2f}",
                    f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}",
                    f"{taxa_sucesso:.1f}"
                ]
            
            comparative_data.append(row)
        
        # Ordenar: algoritmos com erro no final, outros por melhor distância
        def sort_key(row):
            if row[3] == "ERRO":
                return (float('inf'), float(row[1]))  # Erro vai para o final
            return (float(row[3]), float(row[1]))    # Melhor distância, depois tempo
        
        comparative_data.sort(key=sort_key)
        
        table = tabulate(comparative_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        
        # Adicionar ranking
        output.append("\n🥇 RANKING POR PERFORMANCE:")
        output.append("-" * 40)
        for i, row in enumerate(comparative_data, 1):
            if row[3] == "ERRO":
                medal = "❌"
                info = f"Falha em execução | Tempo: {row[1]}s"
            else:
                medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else f"{i}°"
                info = f"Distância: {row[3]} | Tempo: {row[1]}s | Sucesso: {row[6]}%"
            output.append(f"{medal} {row[0]} - {info}")
        
        return "\n".join(output)
    
    def save_detailed_report(self, filename: str = "relatorio_detalhado.txt"):
        """Salva relatório detalhado em arquivo"""
        results_dir = Path(__file__).parent.parent / "results"
        results_dir.mkdir(exist_ok=True)
        filepath = results_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(self.format_detailed_results())


