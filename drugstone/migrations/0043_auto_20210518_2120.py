# Generated by Django 3.0.5 on 2021-05-18 19:20

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0042_auto_20200531_1519'),
    ]

    operations = [
        migrations.CreateModel(
            name='PDIDataset',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='', max_length=128, unique=True)),
                ('link', models.CharField(default='', max_length=128)),
                ('version', models.CharField(default='', max_length=128)),
            ],
        ),
        migrations.CreateModel(
            name='PPIDataset',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='', max_length=128)),
                ('link', models.CharField(default='', max_length=128)),
                ('version', models.CharField(default='', max_length=128)),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='datasetvirus',
            unique_together=None,
        ),
        migrations.AlterUniqueTogether(
            name='proteinviralinteraction',
            unique_together=None,
        ),
        migrations.RemoveField(
            model_name='proteinviralinteraction',
            name='effect',
        ),
        migrations.RemoveField(
            model_name='proteinviralinteraction',
            name='protein',
        ),
        migrations.AlterUniqueTogether(
            name='viralprotein',
            unique_together=None,
        ),
        migrations.RemoveField(
            model_name='viralprotein',
            name='dataset_virus',
        ),
        migrations.RemoveField(
            model_name='drug',
            name='in_literature',
        ),
        migrations.RemoveField(
            model_name='drug',
            name='in_trial',
        ),
        migrations.RemoveField(
            model_name='protein',
            name='closest_effects',
        ),
        migrations.RemoveField(
            model_name='protein',
            name='effects',
        ),
        migrations.AddField(
            model_name='protein',
            name='entrez',
            field=models.CharField(default='', max_length=128, unique=True),
        ),
        migrations.AlterField(
            model_name='tissue',
            name='name',
            field=models.CharField(default='', max_length=128, unique=True),
        ),
        migrations.DeleteModel(
            name='ClosestViralProtein',
        ),
        migrations.DeleteModel(
            name='DatasetVirus',
        ),
        migrations.DeleteModel(
            name='ProteinViralInteraction',
        ),
        migrations.DeleteModel(
            name='ViralProtein',
        ),
        migrations.AddField(
            model_name='proteindruginteraction',
            name='pdi_dataset',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='pdi_dataset_relation', to='drugstone.PDIDataset'),
        ),
        migrations.AddField(
            model_name='proteinproteininteraction',
            name='ppi_dataset',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='ppi_dataset_relation', to='drugstone.PPIDataset'),
        ),
    ]
