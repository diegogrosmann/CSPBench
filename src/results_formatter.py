import statistics
from typing import List, Dict, Any
from tabulate import tabulate
from pathlib import Path

"""
Formatação e salvamento de resultados de execuções de algoritmos CSP.

Classes:
    ResultsFormatter: Armazena, formata e salva resultados detalhados e comparativos.

Métodos:
    add_algorithm_results(...): Adiciona resultados de um algoritmo.
    format_detailed_results(): Formata resultados detalhados.
    save_detailed_report(...): Salva relatório em arquivo.
"""

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
        
        # Extrair informações das bases e algoritmos
        comparative_data = []
        base_algorithms = {}  # Para armazenar algoritmos agrupados por base
        
        for algorithm_full_name, executions in self.results.items():
            # Inicializar variáveis para evitar erros "unbound"
            distancias = []
            tempos = []
            taxa_sucesso = 0.0
            
            # Extrair informação da base e nome do algoritmo
            # Exemplo: "Sintético Pequeno_Base1_CSC" -> Base: "Sintético Pequeno (Base1)", Algoritmo: "CSC"
            parts = algorithm_full_name.split('_')
            
            if len(parts) >= 3 and 'Base' in parts[-2]:
                # Formato: Nome_Base#_Algoritmo
                base_name = '_'.join(parts[:-2])
                base_num = parts[-2]
                algo_name = parts[-1]
                base_display = f"{base_name} ({base_num})"
            else:
                # Formato padrão ou desconhecido
                base_display = "N/A"
                algo_name = algorithm_full_name
            
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
                    base_display,
                    algo_name,
                    f"{statistics.mean(tempos) if tempos else 0:.4f}",
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
                taxa_sucesso = len(valid_executions) / len(executions) * 100 if len(executions) > 0 else 0
                
                row = [
                    base_display,
                    algo_name,
                    f"{statistics.mean(tempos) if tempos else 0:.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    f"{min(distancias) if distancias else 'ERRO'}",
                    f"{statistics.mean(distancias):.2f}" if distancias else "ERRO",
                    f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}" if distancias else "ERRO",
                    f"{taxa_sucesso:.1f}"
                ]
            
            comparative_data.append(row)
            
            # Agrupar por base para o ranking
            if base_display not in base_algorithms:
                base_algorithms[base_display] = []
            
            if valid_executions:
                base_algorithms[base_display].append({
                    'name': algo_name,
                    'dist': min(distancias) if distancias else float('inf'),  # Verificação de segurança
                    'time': statistics.mean(tempos) if tempos else 0,
                    'success_rate': taxa_sucesso if isinstance(taxa_sucesso, (int, float)) else 0  # Garantir que é número
                })
            else:
                base_algorithms[base_display].append({
                    'name': algo_name,
                    'dist': float('inf'),
                    'time': statistics.mean(tempos) if tempos else 0,
                    'success_rate': 0
                })
        
        # Cabeçalhos com base e algoritmo como colunas separadas
        headers = [
            "Base",
            "Algoritmo",
            "Tempo Médio (s)",
            "Desvio Tempo",
            "Melhor Distância",
            "Distância Média", 
            "Desvio Distancia",
            "Taxa Sucesso (%)"
        ]
        
        # Ordenação da tabela: primeiro por base, depois por melhor distância
        def sort_key(row):
            base_name = row[0]
            if row[4] == "ERRO":
                dist = float('inf')
            else:
                dist = float(row[4])
            return (base_name, dist, float(row[2]))  # Base, distância, tempo
        
        comparative_data.sort(key=sort_key)
        
        table = tabulate(comparative_data, headers=headers, tablefmt="grid", stralign="center")
        output.append(table)
        
        # Adicionar ranking POR BASE
        output.append("\n🥇 RANKING POR BASE:")
        output.append("-" * 60)
        
        for base_name, algorithms in sorted(base_algorithms.items()):
            output.append(f"\n💠 BASE: {base_name}")
            output.append("-" * 40)
            
            # Ordenar algoritmos para esta base
            sorted_algorithms = sorted(algorithms, key=lambda x: (x['dist'], x['time']))
            
            # Mostrar ranking por base
            for i, alg in enumerate(sorted_algorithms, 1):
                if alg['dist'] == float('inf'):
                    medal = "❌"
                    info = f"Falha na execução | Tempo: {alg['time']:.4f}s"
                else:
                    medal = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else f"{i}°"
                    info = f"Distância: {alg['dist']} | Tempo: {alg['time']:.4f}s | Sucesso: {alg['success_rate']:.1f}%"
                output.append(f"{medal} {alg['name']} - {info}")
        
        return "\n".join(output)
    
    def save_detailed_report(self, filename: str):
        """
        Salva o relatório detalhado em um arquivo.
        
        Args:
            filename: Nome do arquivo ou caminho completo
        """
        file_path = Path(filename)
        # Garantir que o diretório existe
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(self.format_detailed_results())


